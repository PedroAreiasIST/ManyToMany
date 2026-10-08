// ParallelCutter.cs — concurrent cutting of independent groups of marked edges.
//
// Two groups are independent when no node lies in the elements both change: the
// rows of the split parents, the endpoints and the new nodes. Every write of an edit
// (transpose rows, positions behind an insertion, pattern rows, destinations, facet
// entries of sharers, colours) is then addressed by a node or an element of one group
// only, and every read of an edit is addressed by its own nodes, so the groups can be
// applied by different threads without locks on the mesh structures. The cut runs in
// three phases:
//   plan    (parallel, read-only): shells, children rows, point copies and touched
//           nodes of each group; groups that share a node are merged (union-find) and
//           applied by one thread, which keeps the result exact;
//   reserve (sequential): new node labels in the order of the batch cutter (first
//           appearance in the concatenated groups), one contiguous range of element
//           indices per group, pattern slots bounded from the children rows, private
//           directories (IncrementalMesh.BeginConcurrent);
//   apply   (parallel): each group is cut as by Cut, its children going to its range.
// Node labels, and with them every template decision, are those of the sequential cut
// of the concatenated groups; only the indices of the appended children differ, so the
// rows equal those of the sequential cut as a multiset and every derived structure is
// the one of the resulting mesh. Pattern segment placement and the order of the
// boundary list depend on thread timing; slots and destinations do not.

using System.Diagnostics;

namespace Step2;

public sealed class ParallelCutResult
{
    /// <summary>One record per applied group (merged groups give one record), revertible by Uncut.</summary>
    public required CutRecord[] Groups;

    /// <summary>Input groups merged into another group because they shared a node.</summary>
    public int Merged;

    public double PlanMs, ReserveMs, ApplyMs;
    public long ReservedSlots;
    public int Threads;
}

public static partial class Cutter
{
    private sealed class GroupPlan
    {
        public int Index;
        public List<(int A, int B, int Mid)> Edges = new();
        public List<(int Source, int Node)> Points = new();
        public int[] Parents = [];
        public int[][] Rows = [];
        public List<int[]>[] Children = [];

        /// <summary>Touched nodes (the keys; independence of groups is decided on them).</summary>
        public Dictionary<int, int> Nodes = new();
        public int NewElements, First;
        public long Slots;
    }

    /// <summary>
    ///     Cuts several groups of marked edges concurrently. The result equals
    ///     <see cref="Cut"/> of the concatenated groups up to the indices of the appended
    ///     children. <see cref="ICutObserver.NodeCreated"/> is called sequentially, in label
    ///     order; <see cref="ICutObserver.ElementSplit"/> is called from the writer threads,
    ///     for disjoint children. On an exception the working version is inconsistent and
    ///     must be rolled back to a published version.
    /// </summary>
    public static ParallelCutResult CutParallel(IncrementalMesh m, IReadOnlyList<IReadOnlyList<(int, int)>> groups,
        ICutObserver? obs = null, int threads = 0)
    {
        if (threads <= 0) threads = Environment.ProcessorCount;
        var sw = Stopwatch.StartNew();
        var first = m.NodeCount;
        var e0 = m.ElementCount;

        // ---- labels (sequential): order of first appearance in the concatenated groups
        var mid = new Dictionary<long, int>(EdgeKeyComparer.Instance);
        var plans = new GroupPlan[groups.Count];
        for (var gi = 0; gi < groups.Count; gi++)
        {
            plans[gi] = new GroupPlan { Index = gi };
            foreach (var (p, q) in groups[gi])
            {
                if (p == q) throw new ArgumentException($"degenerate pair ({p},{p})");
                int a = Math.Min(p, q), b = Math.Max(p, q);
                var key = Key(a, b);
                if (mid.ContainsKey(key)) continue;
                var c = m.AddNode();
                mid[key] = c;
                plans[gi].Edges.Add((a, b, c));
                obs?.NodeCreated(c, a, b);
            }
        }

        // ---- plan (parallel, read-only)
        var po = new ParallelOptions { MaxDegreeOfParallelism = threads };
        Parallel.For(0, plans.Length, po, gi => Plan(m, plans[gi], mid));

        // ---- independence: groups sharing a node are merged
        var parent = new int[plans.Length];
        for (var i = 0; i < parent.Length; i++) parent[i] = i;
        int Find(int i)
        {
            while (parent[i] != i) i = parent[i] = parent[parent[i]];
            return i;
        }

        var owner = new Dictionary<int, int>();
        foreach (var pl in plans)
        foreach (var v in pl.Nodes.Keys)
        {
            if (!owner.TryAdd(v, pl.Index))
            {
                int x = Find(owner[v]), y = Find(pl.Index);
                if (x != y) parent[Math.Max(x, y)] = Math.Min(x, y);
            }
        }

        var merged = 0;
        var final = new List<GroupPlan>();
        var members = new Dictionary<int, List<int>>();
        for (var gi = 0; gi < plans.Length; gi++)
        {
            var r = Find(gi);
            if (!members.TryGetValue(r, out var list)) members[r] = list = new List<int>();
            list.Add(gi);
        }

        var replan = new List<GroupPlan>();
        foreach (var (r, list) in members.OrderBy(kv => kv.Key))
        {
            if (list.Count == 1)
            {
                final.Add(plans[r]);
                continue;
            }

            merged += list.Count - 1;
            var g = new GroupPlan { Index = r };
            foreach (var gi in list) g.Edges.AddRange(plans[gi].Edges); // ascending labels
            final.Add(g);
            replan.Add(g);
        }

        Parallel.For(0, replan.Count, po, i => Plan(m, replan[i], mid));
        var planMs = sw.Elapsed.TotalMilliseconds;

        // ---- reserve (sequential)
        var newElements = 0;
        long slots = 0;
        foreach (var g in final)
        {
            g.First = e0 + newElements;
            newElements += g.NewElements;
            slots += g.Slots;
        }

        m.ReserveElements(newElements);
        // a reservation beyond the int range of the pattern arrays: the groups are applied one after another
        // without concurrent mode (same result, arrays grow as needed)
        var concurrent = slots <= int.MaxValue - m.TotalSlots - 2 * 1024;
        if (concurrent) m.BeginConcurrent((int)slots);
        var reserveMs = sw.Elapsed.TotalMilliseconds - planMs;

        // ---- apply (parallel), largest groups first
        var records = new CutRecord[final.Count];
        if (!concurrent)
        {
            for (var i = 0; i < final.Count; i++) records[i] = Apply(m, final[i], obs, first, e0);
            return new ParallelCutResult
            {
                Groups = records, Merged = merged, PlanMs = planMs, ReserveMs = reserveMs,
                ApplyMs = sw.Elapsed.TotalMilliseconds - planMs - reserveMs, ReservedSlots = 0, Threads = 1
            };
        }

        var order = Enumerable.Range(0, final.Count).OrderByDescending(i => final[i].Parents.Length).ToArray();
        var nextGroup = -1;
        try
        {
            var workers = Math.Min(threads, Math.Max(1, final.Count));
            var tasks = new Task[workers];
            for (var w = 0; w < workers; w++)
                tasks[w] = Task.Factory.StartNew(() =>
                {
                    int i;
                    while ((i = Interlocked.Increment(ref nextGroup)) < order.Length)
                        records[order[i]] = Apply(m, final[order[i]], obs, first, e0);
                }, CancellationToken.None, TaskCreationOptions.LongRunning, TaskScheduler.Default);
            Task.WaitAll(tasks);
        }
        finally
        {
            m.EndConcurrent();
        }

        return new ParallelCutResult
        {
            Groups = records, Merged = merged, PlanMs = planMs, ReserveMs = reserveMs,
            ApplyMs = sw.Elapsed.TotalMilliseconds - planMs - reserveMs, ReservedSlots = slots, Threads = threads
        };
    }

    /// <summary>Shells, children, point copies, touched nodes and reservation bounds of one group.</summary>
    private static void Plan(IncrementalMesh m, GroupPlan g, Dictionary<long, int> mid)
    {
        g.Points.Clear();
        g.Nodes.Clear();
        var affected = new HashSet<int>();
        var shell = new List<int>();
        foreach (var (a, b, c) in g.Edges)
        {
            g.Nodes.TryAdd(a, 0);
            g.Nodes.TryAdd(b, 0);
            g.Nodes.TryAdd(c, 0);
            foreach (var v in (ReadOnlySpan<int>)[a, b])
            {
                var tr = m.TransposePairs(v);
                for (var q = 0; q < tr.Length; q += 2)
                    if (m.Row(tr[q]).Length == 1)
                        g.Points.Add((tr[q], c));
            }

            shell.Clear();
            Shell(m, a, b, shell);
            foreach (var e in shell)
                if (m.Row(e).Length >= 2)
                    affected.Add(e);
        }

        g.Parents = affected.ToArray();
        Array.Sort(g.Parents);
        g.Rows = new int[g.Parents.Length][];
        g.Children = new List<int[]>[g.Parents.Length];
        Span<int> mid6 = stackalloc int[6];
        var columns = new Dictionary<int, List<int>>(); // node -> columns of the new rows through it
        void AddColumns(int v, int[] r)
        {
            if (!columns.TryGetValue(v, out var l)) columns[v] = l = new List<int>(16);
            l.AddRange(r);
        }

        var newElements = g.Points.Count;
        for (var i = 0; i < g.Parents.Length; i++)
        {
            var row = m.Row(g.Parents[i]).ToArray();
            var ch = new List<int[]>(8);
            Children(row, mid, mid6, ch);
            g.Rows[i] = row;
            g.Children[i] = ch;
            newElements += ch.Count - 1;
            foreach (var r in ch)
            foreach (var v in r)
                AddColumns(v, r);
        }

        foreach (var (_, c) in g.Points) AddColumns(c, [c]);

        // a row ends with at most its current columns and the distinct new columns not among them
        long slots = 0;
        var distinct = new HashSet<int>();
        foreach (var (v, list) in columns)
        {
            g.Nodes.TryAdd(v, 0);
            distinct.Clear();
            var cur = m.PatternColumns(v);
            foreach (var c in list)
                if (cur.IndexOf(c) < 0)
                    distinct.Add(c);
            slots += m.SlotsToGrow(v, cur.Length + distinct.Count);
        }

        g.Slots = slots;
        g.NewElements = newElements;
    }

    private static CutRecord Apply(IncrementalMesh m, GroupPlan g, ICutObserver? obs, int first, int e0)
    {
        var next = g.First;
        var points = new (int, int)[g.Points.Count];
        for (var i = 0; i < g.Points.Count; i++)
        {
            var (src, c) = g.Points[i];
            m.CreateElement(next, [c]);
            obs?.ElementSplit(src, next);
            points[i] = (next++, src);
        }

        var splits = new (int, int[], int[])[g.Parents.Length];
        for (var i = 0; i < g.Parents.Length; i++)
        {
            var e = g.Parents[i];
            var ch = g.Children[i];
            var idx = new int[ch.Count];
            idx[0] = e;
            m.ReplaceElement(e, ch[0]);
            obs?.ElementSplit(e, e);
            for (var t = 1; t < ch.Count; t++)
            {
                idx[t] = next++;
                m.CreateElement(idx[t], ch[t]);
                obs?.ElementSplit(e, idx[t]);
            }

            splits[i] = (e, g.Rows[i], idx);
        }

        if (next != g.First + g.NewElements) throw new InvalidOperationException("element reservation mismatch");
        return new CutRecord
        {
            Edges = g.Edges.ToArray(), Splits = splits, Points = points, FirstNewNode = first, ElementsBefore = e0,
            ShellElements = g.Parents.Length
        };
    }
}
