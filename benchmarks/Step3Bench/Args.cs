// Args.cs — key=value options of the driver modes.
//
// Copyright (C) 2026 Pedro Miguel de Almeida Areias
//
// This file is part of ManyToMany.
//
// ManyToMany is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ManyToMany is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ManyToMany. If not, see <https://www.gnu.org/licenses/>.

using System.Globalization;

namespace Step3.Bench;

public sealed class Args
{
    private readonly Dictionary<string, string> _kv = new(StringComparer.OrdinalIgnoreCase);
    private readonly HashSet<string> _flags = new(StringComparer.OrdinalIgnoreCase);
    public readonly List<string> Positional = new();

    public Args(IEnumerable<string> a)
    {
        foreach (var s in a)
        {
            var i = s.IndexOf('=');
            if (i > 0) _kv[s[..i]] = s[(i + 1)..];
            else if (double.TryParse(s, NumberStyles.Float, CultureInfo.InvariantCulture, out _)) Positional.Add(s);
            else _flags.Add(s.TrimStart('-'));
        }
    }

    public double D(string key, double def) => _kv.TryGetValue(key, out var v) ? double.Parse(v, CultureInfo.InvariantCulture) : def;
    public int I(string key, int def) => _kv.TryGetValue(key, out var v) ? int.Parse(v, CultureInfo.InvariantCulture) : def;
    public string S(string key, string def) => _kv.TryGetValue(key, out var v) ? v : def;
    public bool Flag(string key) => _flags.Contains(key) || (_kv.TryGetValue(key, out var v) && (v == "1" || v.Equals("true", StringComparison.OrdinalIgnoreCase)));
    public bool Has(string key) => _kv.ContainsKey(key);
}
