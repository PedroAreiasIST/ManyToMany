# Topology Library: Public API Reference and User Guide

**Author: Pedro Areias**
**Public API Documentation for Computational Mechanics Applications**

---

## Preface

This document provides complete documentation for the **public API** of the Topology library. Every class, struct, enum, method, and property described here is part of the public contract available to users. Internal implementation details are omitted.

**Prerequisites:** Familiarity with C# programming and basic understanding of finite element analysis and mesh data structures.

**Requirements:** Modern C# features including collection expressions, primary constructors, and span-based APIs are utilized throughout.

---

## Table of Contents

### Part I: Getting Started

1. [Introduction](#1-introduction)
2. [Core Types](#2-core-types)
3. [Creating a Topology](#3-creating-a-topology)
4. [Entity Operations](#4-entity-operations)
5. [Data Management](#5-data-management)

### Part II: Queries and Connectivity

6. [Connectivity Queries](#6-connectivity-queries)
7. [Neighbor and Adjacency Queries](#7-neighbor-and-adjacency-queries)
8. [Set Operations](#8-set-operations)
9. [Multi-Type Connectivity](#9-multi-type-connectivity)

### Part III: Advanced Features

10. [Sub-Entity Extraction](#10-sub-entity-extraction)
11. [Boundary Detection](#11-boundary-detection)
12. [Graph Algorithms](#12-graph-algorithms)
13. [Symmetry and Canonical Forms](#13-symmetry-and-canonical-forms)
14. [Graph Coloring](#14-graph-coloring)
15. [Dual Graph Construction](#15-dual-graph-construction)
16. [Bandwidth Reduction](#16-bandwidth-reduction)
17. [Sparse Matrix Patterns](#17-sparse-matrix-patterns)

### Part IV: Smart Handles and Traversal

18. [Smart Entity Handles](#18-smart-entity-handles)
19. [Traversal and Components](#19-traversal-and-components)

### Part V: Operations and Performance

20. [Batch Operations](#20-batch-operations)
21. [Compression and Optimization](#21-compression-and-optimization)
22. [Mesh Merging and Extraction](#22-mesh-merging-and-extraction)
23. [Serialization](#23-serialization)
24. [Performance Configuration](#24-performance-configuration)

### Part VI: Supporting Types

25. [O2M — One-to-Many Sparse Structure](#25-o2m--one-to-many-sparse-structure)
26. [M2M — Thread-Safe Many-to-Many](#26-m2m--thread-safe-many-to-many)
27. [MM2M — Multi-Type Manager](#27-mm2m--multi-type-manager)
28. [Utility Functions](#28-utility-functions)

### Appendices

A. [API Quick Reference](#appendix-a-api-quick-reference)
B. [Performance Characteristics](#appendix-b-performance-characteristics)
C. [Common Patterns](#appendix-c-common-patterns)

---

# Part I: Getting Started

## 1. Introduction

The Topology library provides a type-safe, high-performance mesh data structure for finite element analysis and computational mechanics. It manages connectivity relationships between entities (nodes, edges, faces, elements) with compile-time type checking and zero-allocation query paths.

**Key capabilities:**

- Compile-time type safety via generic type maps
- Zero-allocation span-based APIs for hot paths
- Thread-safe operations with reader-writer locking
- Symmetry-aware deduplication
- Graph algorithms (BFS, Dijkstra, connected components)
- Sparse matrix pattern extraction for FEA assembly

```csharp
using Numerical;

// Create a mesh with nodes and elements
using var mesh = Topology.New<Node, Element>();

// Add nodes
int n0 = mesh.Add<Node, Point>(new Point(0, 0, 0));
int n1 = mesh.Add<Node, Point>(new Point(1, 0, 0));
int n2 = mesh.Add<Node, Point>(new Point(0, 1, 0));

// Add a triangular element
int elem = mesh.Add<Element, Node>(n0, n1, n2);

// Query connectivity
IReadOnlyList<int> nodes = mesh.NodesOf<Element, Node>(elem);
IReadOnlyList<int> elements = mesh.ElementsAt<Element, Node>(n0);
```

---

## 2. Core Types

### 2.1 ITypeMap Interface

Defines the type-to-index mapping used by `Topology<TTypes>`.

```csharp
public interface ITypeMap
{
    int Count { get; }
    int IndexOf<T>();
    bool TryIndexOf<T>(out int index);
}
```

### 2.2 TypeMap Classes

Pre-built `ITypeMap` implementations for 2 to 25 types. Each is a sealed class.

```csharp
public sealed class TypeMap<T0, T1> : ITypeMap { ... }
public sealed class TypeMap<T0, T1, T2> : ITypeMap { ... }
public sealed class TypeMap<T0, T1, T2, T3> : ITypeMap { ... }
// ... up to TypeMap<T0, ..., T24>
```

Type indices are assigned in declaration order:

```csharp
// In TypeMap<Node, Edge, Face, Element>:
// Node  → index 0
// Edge  → index 1
// Face  → index 2
// Element → index 3
```

### 2.3 ResultOrder Enum

```csharp
public enum ResultOrder
{
    Unordered = 0,  // Insertion order (fastest)
    Sorted = 1      // Sorted by (type index, entity index)
}
```

### 2.4 SubEntityDefinition Struct

Defines local node index combinations for extracting sub-entities (edges, faces) from parent elements.

```csharp
public readonly struct SubEntityDefinition
{
    public readonly int[][] LocalNodeIndices;

    public SubEntityDefinition(int[][] localNodeIndices);
    public static SubEntityDefinition FromEdges(params (int, int)[] edges);
    public static SubEntityDefinition FromFaces(params (int, int, int)[] faces);
    public static SubEntityDefinition FromQuadFaces(params (int, int, int, int)[] faces);
}
```

**Example — triangle edge definition:**

```csharp
var triEdges = SubEntityDefinition.FromEdges(
    (0, 1), (1, 2), (2, 0)
);

var tetFaces = SubEntityDefinition.FromFaces(
    (0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)
);
```

---

## 3. Creating a Topology

### 3.1 Static Factory: Topology Class

The non-generic `Topology` class provides factory methods for creating typed topologies:

```csharp
public static class Topology
{
    public static Topology<TypeMap<T0, T1>> New<T0, T1>();
    public static Topology<TypeMap<T0, T1, T2>> New<T0, T1, T2>();
    public static Topology<TypeMap<T0, T1, T2, T3>> New<T0, T1, T2, T3>();
    // ... up to 25 type parameters
}
```

**Usage:**

```csharp
// 2D mesh with nodes and elements
using var mesh = Topology.New<Node, Element>();

// Full 3D mesh hierarchy
using var mesh3D = Topology.New<Node, Edge, Face, Element>();

// Mixed element types
using var mixed = Topology.New<Node, Edge, Face, Tet, Hex, Prism, Pyramid>();
```

### 3.2 Topology\<TTypes\> Constructor

```csharp
public class Topology<TTypes> : IDisposable where TTypes : ITypeMap, new()
{
    public Topology();
}
```

Directly instantiable when using a custom `ITypeMap`:

```csharp
using var mesh = new Topology<TypeMap<Node, Element>>();
```

### 3.3 Disposal

`Topology<TTypes>` implements `IDisposable`. Always use `using` statements:

```csharp
using var mesh = Topology.New<Node, Element>();
// mesh is disposed at end of scope
```

### 3.4 Cloning and Read-Only Views

```csharp
public Topology<TTypes> Clone();
public ReadOnlyTopology<TTypes> AsReadOnly();
```

```csharp
// Deep copy
var copy = mesh.Clone();

// Read-only view (query-only, no modification methods)
ReadOnlyTopology<TTypes> readOnly = mesh.AsReadOnly();
```

---

## 4. Entity Operations

### 4.1 Adding Standalone Entities

```csharp
public int Add<TNode>();
public int Add<TNode, TData>(TData data);
```

```csharp
int n0 = mesh.Add<Node>();                              // No data
int n1 = mesh.Add<Node, Point>(new Point(1, 0, 0));     // With data
```

### 4.2 Adding Connected Entities

```csharp
public int Add<TElement, TNode>(params int[] connectedNodes);
public int Add<TElement, TNode>(ReadOnlySpan<int> connectedNodes);
public int Add<TElement, TNode, TData>(TData data, params int[] connectedNodes);
```

```csharp
int elem = mesh.Add<Element, Node>(n0, n1, n2);

var material = new Material(E: 210e9, Nu: 0.3);
int elem = mesh.Add<Element, Node, Material>(material, n0, n1, n2, n3);

// Span overload (zero allocation)
ReadOnlySpan<int> nodes = [n0, n1, n2];
int elem = mesh.Add<Element, Node>(nodes);
```

### 4.3 Unique Addition (Deduplication)

Requires symmetry configuration. Returns existing index if equivalent element exists.

```csharp
public (int Index, bool WasNew) AddUnique<TElement, TNode>(params int[] connectedNodes);
public (int Index, bool WasNew) AddUnique<TElement, TNode>(ReadOnlySpan<int> connectedNodes);
public (int Index, bool WasNew) AddUnique<TElement, TNode, TData>(TData data, params int[] connectedNodes);
```

```csharp
mesh.WithSymmetry<Edge>(Symmetry.Dihedral(2));

var (idx1, isNew1) = mesh.AddUnique<Edge, Node>(n0, n1);  // Creates edge, isNew=true
var (idx2, isNew2) = mesh.AddUnique<Edge, Node>(n1, n0);  // Returns same, isNew=false
// idx1 == idx2
```

### 4.4 Bulk Addition

```csharp
public int[] AddRange<TNode, TData>(IEnumerable<TData> dataItems);
public int[] AddRange<TNode, TData>(ReadOnlySpan<TData> dataItems);
public int[] AddRange<TElement, TNode>(IEnumerable<int[]> connectivityList);
public (int Index, bool WasNew)[] AddRangeUnique<TElement, TNode>(IEnumerable<int[]> connectivityList);
public int[] AddRangeParallel<TElement, TNode>(int[][] connectivityList, int minParallelCount = 10000);
```

```csharp
// Add multiple nodes with data
Point[] coords = [new(0,0,0), new(1,0,0), new(0,1,0)];
int[] nodeIds = mesh.AddRange<Node, Point>(coords);

// Add multiple elements
int[][] connectivity = [[0,1,2,3], [1,4,5,2], [2,5,6,3]];
int[] elemIds = mesh.AddRange<Element, Node>(connectivity);

// Parallel addition (large meshes)
int[] elemIds = mesh.AddRangeParallel<Element, Node>(connectivity, minParallelCount: 5000);
```

### 4.5 Existence and Lookup

```csharp
public bool Exists<TElement>(params int[] nodes);
public int Find<TElement>(params int[] nodes);
```

```csharp
bool exists = mesh.Exists<Edge>(n0, n1);      // true/false
int edgeIdx = mesh.Find<Edge>(n0, n1);         // Index or -1
```

### 4.6 Removal

Entities are marked for removal and physically removed during `Compress()`.

```csharp
public void Remove<TEntity>(int index);
public void RemoveRange<TEntity>(IEnumerable<int> indices);
```

```csharp
mesh.Remove<Node>(5);
mesh.RemoveRange<Node>([5, 10, 15]);
mesh.Compress();  // Apply pending removals
```

### 4.7 In-Place Element Modification

```csharp
public void AddNodeToElement<TElement, TNode>(int element, int node);
public bool RemoveNodeFromElement<TElement, TNode>(int element, int node);
public void ReplaceElementNodes<TElement, TNode>(int element, params int[] newNodes);
public void ClearElement<TElement, TNode>(int element);
```

```csharp
// Add node to existing element
mesh.AddNodeToElement<Element, Node>(elemIdx, newNodeIdx);

// Replace all nodes (e.g. during mesh refinement)
mesh.ReplaceElementNodes<Element, Node>(elemIdx, n5, n6, n7, n8);

// Remove specific node
bool removed = mesh.RemoveNodeFromElement<Element, Node>(elemIdx, oldNodeIdx);
```

### 4.8 Counting

```csharp
public int Count<TEntity>();
public int CountActive<TEntity>();
public int CountRelated<TEntity, TRelated>(int entityIndex);
public int CountIncident<TElement, TNode>(int nodeIndex);
```

```csharp
int nodeCount = mesh.Count<Node>();
int activeNodes = mesh.CountActive<Node>();       // Excludes marked-for-removal
int nodesInElem = mesh.CountRelated<Element, Node>(elemIdx);
int elemsAtNode = mesh.CountIncident<Element, Node>(nodeIdx);
```

### 4.9 Active and Marked Entities

```csharp
public List<int> GetActive<TEntity>();
public IReadOnlySet<int> GetMarkedForRemoval<TEntity>();
```

---

## 5. Data Management

### 5.1 Get and Set

```csharp
public TData Get<TEntity, TData>(int index);
public bool TryGet<TEntity, TData>(int index, out TData? value);
public void Set<TEntity, TData>(int index, TData value);
public void SetRange<TEntity, TData>(int startIndex, ReadOnlySpan<TData> values);
```

```csharp
mesh.Set<Node, Point>(nodeIdx, new Point(1.0, 2.0, 3.0));
Point coord = mesh.Get<Node, Point>(nodeIdx);

if (mesh.TryGet<Node, Point>(nodeIdx, out Point p))
    Console.WriteLine($"Node at {p}");

// Bulk set
Point[] coords = GenerateCoordinates(100);
mesh.SetRange<Node, Point>(0, coords);
```

### 5.2 Retrieval and Iteration

```csharp
public IReadOnlyList<TData> All<TEntity, TData>();
public IEnumerable<(int Index, TData Data)> Each<TEntity, TData>();
public IEnumerable<int> Each<TEntity>();
public void ForEach<TEntity, TData>(Action<int, TData> action);
public void ParallelForEach<TEntity, TData>(Action<int, TData> action, int minParallelCount = 1000);
```

```csharp
// Iterate indices
foreach (int idx in mesh.Each<Node>())
    Console.WriteLine($"Node {idx}");

// Iterate with data
foreach (var (idx, coord) in mesh.Each<Node, Point>())
    Console.WriteLine($"Node {idx}: {coord}");

// Callback-based (avoids iterator allocation)
mesh.ForEach<Node, Point>((idx, coord) => ProcessNode(idx, coord));

// Parallel iteration
mesh.ParallelForEach<Element, Material>((idx, mat) =>
{
    var Ke = ComputeElementStiffness(idx, mat);
    StoreLocal(idx, Ke);
}, minParallelCount: 1000);
```

---

# Part II: Queries and Connectivity

## 6. Connectivity Queries

### 6.1 Forward Queries (Element → Nodes)

```csharp
public IReadOnlyList<int> NodesOf<TElement, TNode>(int element);
public void WithNodesOf<TElement, TNode>(int element, Action<ReadOnlySpan<int>> action);
public TResult WithNodesOf<TElement, TNode, TResult>(int element, Func<ReadOnlySpan<int>, TResult> func);
public void WithNodesSpan<TElement, TNode>(int element, Action<ReadOnlySpan<int>> action);
public TResult WithNodesSpan<TElement, TNode, TResult>(int element, Func<ReadOnlySpan<int>, TResult> func);
```

```csharp
// List-based (allocates defensive copy)
IReadOnlyList<int> nodes = mesh.NodesOf<Element, Node>(elemIdx);

// Zero-copy span access (preferred for hot paths)
mesh.WithNodesSpan<Element, Node>(elemIdx, nodes =>
{
    for (int i = 0; i < nodes.Length; i++)
        Process(nodes[i]);
});

// Zero-copy with return value
double area = mesh.WithNodesSpan<Element, Node, double>(elemIdx, nodes =>
{
    return ComputeArea(nodes);
});
```

### 6.2 Reverse Queries (Node → Elements)

The reverse (transpose) is computed on first access and cached.

```csharp
public IReadOnlyList<int> ElementsAt<TElement, TNode>(int node);
```

```csharp
IReadOnlyList<int> elements = mesh.ElementsAt<Element, Node>(nodeIdx);
```

### 6.3 Local Index Lookup

```csharp
public int GetLocalNodeIndex<TElement, TNode>(int element, int node);
public IReadOnlyList<int> GetElementPositions<TElement, TNode>(int element);
public IReadOnlyList<int> GetNodePositions<TElement, TNode>(int node);
```

```csharp
// Which local position does node occupy in element?
int localIdx = mesh.GetLocalNodeIndex<Element, Node>(elemIdx, nodeIdx);
```

---

## 7. Neighbor and Adjacency Queries

### 7.1 Element Neighbors

```csharp
public IReadOnlyList<int> Neighbors<TElement, TNode>(int element, bool sorted = true);
public IEnumerable<int> EnumerateNeighbors<TElement, TNode>(int element);
public List<int> GetDirectNeighbors<TEntity, TRelated>(int entityIndex, bool includeSelf = false, bool sorted = true);
public List<int> GetElementNeighbors<TElement, TNode>(int element, bool sorted = true);
public List<int> GetNodeNeighbors<TElement, TNode>(int node, bool sorted = true);
```

```csharp
// All elements sharing at least one node
var neighbors = mesh.Neighbors<Element, Node>(elemIdx);

// Memory-efficient lazy enumeration
foreach (int neighbor in mesh.EnumerateNeighbors<Element, Node>(elemIdx))
    Process(neighbor);

// Direct neighbors (generic version)
var direct = mesh.GetDirectNeighbors<Element, Node>(elemIdx, sorted: true);
```

### 7.2 Shared-Count Queries

```csharp
public List<int> GetEntitiesWithSharedCount<TEntity, TRelated>(int entityIndex, int exactCount, bool includeSelf = false);
public List<int> GetEntitiesWithMinSharedCount<TEntity, TRelated>(int entityIndex, int minCount = 1, bool includeSelf = false);
public List<(int EntityIndex, int SharedCount)> GetWeightedNeighbors<TEntity, TRelated>(int entityIndex, int minCount = 1, bool includeSelf = false);
```

```csharp
// Elements sharing exactly 2 nodes (face-adjacent in 3D)
var faceNeighbors = mesh.GetEntitiesWithSharedCount<Element, Node>(elemIdx, exactCount: 2);

// Elements sharing at least 1 node
var anyShared = mesh.GetEntitiesWithMinSharedCount<Element, Node>(elemIdx, minCount: 1);

// Neighbors with shared-node counts
var weighted = mesh.GetWeightedNeighbors<Element, Node>(elemIdx);
foreach (var (neighbor, sharedCount) in weighted)
    Console.WriteLine($"Element {neighbor} shares {sharedCount} nodes");
```

### 7.3 K-Hop Neighborhood

```csharp
public Dictionary<int, int> GetKHopNeighborhood<TEntity, TRelated>(
    int seedEntity, int k, int minSharedForConnection = 1, bool includeSeed = true);
public List<int> GetEntitiesAtDistance<TEntity, TRelated>(
    int seedEntity, int k, int minSharedForConnection = 1);
```

```csharp
// All elements within 2 hops (value = distance)
var neighborhood = mesh.GetKHopNeighborhood<Element, Node>(elemIdx, k: 2);

// Elements at exactly distance 3
var ring = mesh.GetEntitiesAtDistance<Element, Node>(elemIdx, k: 3);
```

### 7.4 Element Lookup by Nodes

```csharp
public List<int> GetElementsWithNodes<TElement, TNode>(List<int> nodes);
public List<int> GetElementsContainingAnyNode<TElement, TNode>(List<int> nodes);
public List<int> GetElementsFromNodes<TElement, TNode>(List<int> nodes);
public List<int> ElementsContainingAllNodes<TElement, TNode>(params int[] nodes);
```

```csharp
// Elements with exactly these nodes
var exact = mesh.GetElementsWithNodes<Element, Node>([n0, n1, n2]);

// Elements containing any of the given nodes
var any = mesh.GetElementsContainingAnyNode<Element, Node>([n0, n1]);

// Elements connected to ALL given nodes
var all = mesh.GetElementsFromNodes<Element, Node>([n0, n1]);
```

---

## 8. Set Operations

```csharp
public IReadOnlyList<int> ElementsAtAll<TElement, TNode>(params int[] nodes);
public IReadOnlyList<int> ElementsAtAny<TElement, TNode>(params int[] nodes);
public IReadOnlyList<int> ElementsAtExcluding<TElement, TNode>(int[] include, int[] exclude);
```

```csharp
// Intersection: elements at ALL specified nodes
var shared = mesh.ElementsAtAll<Element, Node>(n0, n1);

// Union: elements at ANY specified node
var any = mesh.ElementsAtAny<Element, Node>(n0, n1, n2);

// Difference: elements at 'include' nodes but NOT at 'exclude' nodes
var diff = mesh.ElementsAtExcluding<Element, Node>([n0, n1], [n2]);
```

---

## 9. Multi-Type Connectivity

### 9.1 Cross-Type Queries

```csharp
public List<(int TypeIndex, int EntityIndex)> MultiTypeDFS<TNode>(int nodeIndex);
public List<(int TypeIndex, int EntityIndex)> GetAllEntitiesAtNode<TNode>(int nodeIndex);
public List<(int TypeIndex, int EntityIndex)> GetAllEntitiesAtNode<TNode>(int nodeIndex, ResultOrder order);
public List<(int TypeIndex, int NodeIndex)> GetAllNodesOfEntity<TEntity>(int entityIndex);
public List<(int TypeIndex, int NodeIndex)> GetAllNodesOfEntity<TEntity>(int entityIndex, ResultOrder order);
```

```csharp
// All entities (edges, faces, elements) connected to a node
var entities = mesh.GetAllEntitiesAtNode<Node>(nodeIdx);
foreach (var (typeIdx, entIdx) in entities)
    Console.WriteLine($"Type {typeIdx}, Entity {entIdx}");

// All nodes of an entity across all node types
var nodes = mesh.GetAllNodesOfEntity<Element>(elemIdx, ResultOrder.Sorted);

// DFS across all types from a node
var reachable = mesh.MultiTypeDFS<Node>(nodeIdx);
```

### 9.2 Type Dependency Analysis

```csharp
public List<int> GetTypeTopologicalOrder();
public bool IsTypeHierarchyAcyclic();
public IReadOnlyList<int> GetTypeDependencyOrder();
public bool AreTypeDependenciesAcyclic();
public IReadOnlyList<int> GetDependencies<TEntity>();
public IReadOnlyList<int> GetDependents<TEntity>();
```

```csharp
var order = mesh.GetTypeDependencyOrder();
bool acyclic = mesh.AreTypeDependenciesAcyclic();

// What types does Element depend on?
var deps = mesh.GetDependencies<Element>();

// What types depend on Node?
var dependents = mesh.GetDependents<Node>();
```

### 9.3 Algebraic Connectivity

```csharp
public O2M ComputeTransitiveConnectivity<TEntity, TRelated>();
public O2M GetDualStructure<TEntity, TRelated>();
public O2M ComputeRelatedToRelatedConnectivity<TEntity, TRelated>();
public bool IsAcyclic<TEntity, TRelated>();
public List<int> GetTopologicalOrder<TEntity, TRelated>();
```

```csharp
// Element-to-element via shared nodes
O2M elemToElem = mesh.ComputeRelatedToRelatedConnectivity<Element, Node>();

// Transitive closure
O2M transitive = mesh.ComputeTransitiveConnectivity<Element, Node>();

// Dual structure (transpose)
O2M dual = mesh.GetDualStructure<Element, Node>();
```

### 9.4 Graph Construction

```csharp
public O2M GetElementToElementGraph<TElement, TNode>();
public O2M GetNodeToNodeGraph<TElement, TNode>();
```

```csharp
O2M elemGraph = mesh.GetElementToElementGraph<Element, Node>();
O2M nodeGraph = mesh.GetNodeToNodeGraph<Element, Node>();
```

---

# Part III: Advanced Features

## 10. Sub-Entity Extraction

```csharp
public (int TotalExtracted, int UniqueAdded, int DuplicatesSkipped)
    DiscoverSubEntities<TElement, TSubEntity, TNode>(SubEntityDefinition definition, bool addUnique = true);
public List<int> ElementsSharingSubEntity<TParent, TSubEntity, TNode>(int subEntityIndex);
public int CountElementsSharingSubEntity<TParent, TSubEntity, TNode>(int subEntityIndex);
```

```csharp
using var mesh = Topology.New<Node, Edge, Face, Element>();
mesh.WithSymmetry<Edge>(Symmetry.Dihedral(2));

// Extract edges from triangular elements
var triEdges = SubEntityDefinition.FromEdges((0,1), (1,2), (2,0));
var (total, unique, dupes) = mesh.DiscoverSubEntities<Element, Edge, Node>(triEdges);
Console.WriteLine($"Extracted {total} edges, {unique} unique, {dupes} duplicates");

// Which elements share a given edge?
var sharing = mesh.ElementsSharingSubEntity<Element, Edge, Node>(edgeIdx);
```

---

## 11. Boundary Detection

```csharp
public HashSet<int> FindBoundaryNodes<TElement, TNode>(int nodesPerBoundaryFacet);
public HashSet<int> FindBoundaryElements<TElement, TNode>(int nodesPerBoundaryFacet);
public List<int[]> ExtractBoundaryFacets<TElement, TNode>(int nodesPerBoundaryFacet);
public List<(int[] Nodes, int Element1, int Element2)> FindInternalFacets<TElement, TNode>(int nodesPerFacet);
```

```csharp
// For triangular meshes, boundary edges have 2 nodes
var boundaryNodes = mesh.FindBoundaryNodes<Element, Node>(nodesPerBoundaryFacet: 2);
var boundaryElems = mesh.FindBoundaryElements<Element, Node>(nodesPerBoundaryFacet: 2);
var boundaryEdges = mesh.ExtractBoundaryFacets<Element, Node>(nodesPerBoundaryFacet: 2);
var internalEdges = mesh.FindInternalFacets<Element, Node>(nodesPerFacet: 2);
```

### 11.1 Sub-Entity Boundary Detection

```csharp
public SubEntityBoundaryResult DetectSubEntityBoundary<TParent, TSubEntity, TNode>();
public List<int> GetBoundarySubEntities<TParent, TSubEntity, TNode>();
public List<int> GetInteriorSubEntities<TParent, TSubEntity, TNode>();
public bool IsSubEntityOnBoundary<TParent, TSubEntity, TNode>(int subEntityIndex);
public List<int> DetectNonManifoldSubEntities<TParent, TSubEntity, TNode>();
```

```csharp
var result = mesh.DetectSubEntityBoundary<Element, Edge, Node>();
Console.WriteLine($"Boundary edges: {result.BoundaryCount}");
Console.WriteLine($"Interior edges: {result.InteriorCount}");

bool isBoundary = mesh.IsSubEntityOnBoundary<Element, Edge, Node>(edgeIdx);
var nonManifold = mesh.DetectNonManifoldSubEntities<Element, Edge, Node>();
```

### SubEntityBoundaryResult

```csharp
public readonly struct SubEntityBoundaryResult
{
    public List<int> BoundaryIndices { get; }
    public List<int> InteriorIndices { get; }
    public int[] IncidenceCounts { get; }
    public int BoundaryCount { get; }
    public int InteriorCount { get; }
}
```

---

## 12. Graph Algorithms

### 12.1 BFS and Dijkstra on Topology

```csharp
public List<int> BreadthFirstSearch<TEntity>(int startEntity, Action<int, int>? visitor = null);
public Dictionary<int, int> BreadthFirstDistances<TEntity>(int startEntity);
public Dictionary<int, (double Distance, int Predecessor)> DijkstraShortestPaths<TEntity>(
    int startEntity, Func<int, int, int, double> edgeWeight);
```

```csharp
// BFS traversal
var visited = mesh.BreadthFirstSearch<Element>(startElem);

// BFS with visitor callback
mesh.BreadthFirstSearch<Element>(startElem, (entity, depth) =>
{
    Console.WriteLine($"Visited element {entity} at depth {depth}");
});

// BFS distances
var distances = mesh.BreadthFirstDistances<Element>(startElem);

// Dijkstra shortest paths
var paths = mesh.DijkstraShortestPaths<Element>(startElem,
    (from, to, edgeIndex) => ComputeEdgeWeight(from, to));
```

### 12.2 Connected Components

```csharp
public List<List<int>> FindConnectedComponents<TEntity, TRelated>(int minSharedForConnection = 1);
public ConnectivityStatistics GetConnectivityStatistics<TEntity, TRelated>(params int[] sharedCounts);
public Dictionary<int, int> ComputeDegrees<TEntity, TRelated>();
```

```csharp
var components = mesh.FindConnectedComponents<Element, Node>();
Console.WriteLine($"Found {components.Count} connected components");

var degrees = mesh.ComputeDegrees<Element, Node>();
var stats = mesh.GetConnectivityStatistics<Element, Node>(1, 2, 3);
```

### ConnectivityStatistics

```csharp
public sealed class ConnectivityStatistics
{
    public int EntityCount { get; }
    public Dictionary<int, List<int>> NeighborCountsByLevel { get; }
    public double GetAverageNeighbors(int sharedCount);
    public int GetMinNeighbors(int sharedCount);
    public int GetMaxNeighbors(int sharedCount);
}
```

### 12.3 Validation and Duplicates

```csharp
public bool ValidateStructure();
public async Task<bool> ValidateStructureAsync(CancellationToken cancellationToken = default);
public ValidationResult ValidateIntegrity<TElement, TNode>();
public List<int> GetDuplicates<TEntity>();
public Dictionary<int, List<int>> GetAllDuplicates();
```

### ValidationResult

```csharp
public readonly struct ValidationResult
{
    public bool IsValid { get; }
    public IReadOnlyList<string> Errors { get; }
}
```

### 12.4 Element Ordering

```csharp
public List<int> GetTopologicalOrder<TEntity>();
public List<int> GetSortOrder<TEntity>();
```

---

## 13. Symmetry and Canonical Forms

### 13.1 Symmetry Class

```csharp
public sealed class Symmetry
{
    public Symmetry(List<List<int>> permutations);

    // Properties
    public int NodeCount { get; }
    public int GroupSize { get; }
    public IReadOnlyList<IReadOnlyList<int>> Permutations { get; }

    // Canonical form computation
    public List<int> Canonical(params int[] nodes);
    public List<int> Canonical(List<int> nodes);
    public List<int> Canonical(IReadOnlyList<int> nodes);
    public void CanonicalSpan(ReadOnlySpan<int> nodes, Span<int> destination);

    // Equivalence checking
    public bool AreEquivalent(int[] a, int[] b);
    public bool AreEquivalent(List<int> a, List<int> b);
    public bool AreEquivalent(IReadOnlyList<int> a, IReadOnlyList<int> b);

    // Permutation application
    public List<int> Apply(IReadOnlyList<int> nodes, int permutationIndex);

    // Factory methods
    public static Symmetry Identity(int nodeCount);
    public static Symmetry Cyclic(int n);
    public static Symmetry Dihedral(int n);
    public static Symmetry Full(int n);
    public static Symmetry FromGenerators(List<List<int>> generators);
}
```

### 13.2 Predefined Symmetries

| Factory | Description | Order | Use Case |
|---------|-------------|-------|----------|
| `Identity(n)` | No symmetry | 1 | Directed edges |
| `Cyclic(n)` | Rotations only | n | Oriented faces |
| `Dihedral(n)` | Rotations + reflections | 2n | Undirected edges (n=2), quad faces (n=4) |
| `Full(n)` | All permutations | n! | Tetrahedra (n=4) |
| `FromGenerators(g)` | Custom group | varies | Custom element types |

### 13.3 Configuration

```csharp
public Topology<TTypes> WithSymmetry<TElement>(Symmetry symmetry);
public Symmetry? GetSymmetry<TElement>();
```

```csharp
// Configure before adding entities
mesh.WithSymmetry<Edge>(Symmetry.Dihedral(2));
mesh.WithSymmetry<Triangle>(Symmetry.Cyclic(3));
mesh.WithSymmetry<Tet>(Symmetry.Full(4));

// Fluent chaining
var mesh = Topology.New<Node, Edge, Face>()
    .WithSymmetry<Edge>(Symmetry.Dihedral(2))
    .WithSymmetry<Face>(Symmetry.Dihedral(3));
```

### 13.4 Usage Example

```csharp
using var mesh = Topology.New<Node, Edge>();
mesh.WithSymmetry<Edge>(Symmetry.Dihedral(2));

int n0 = mesh.Add<Node>(), n1 = mesh.Add<Node>();

var (idx, isNew) = mesh.AddUnique<Edge, Node>(n0, n1);   // Creates
var (idx2, isNew2) = mesh.AddUnique<Edge, Node>(n1, n0); // Returns same
// idx == idx2, isNew2 == false

bool equiv = Symmetry.Dihedral(2).AreEquivalent([0, 1], [1, 0]); // true
```

---

## 14. Graph Coloring

```csharp
public int[] ComputeElementColoring<TElement, TNode>();
public List<List<int>> GetColorGroups<TElement, TNode>();
public ColoringStatistics GetColoringStatistics<TElement, TNode>();
```

```csharp
// Color elements so no two neighbors share a color
int[] colors = mesh.ComputeElementColoring<Element, Node>();

// Get groups by color for parallel assembly
var groups = mesh.GetColorGroups<Element, Node>();
foreach (var group in groups)
{
    // All elements in a group can be assembled in parallel
    Parallel.ForEach(group, elem => AssembleElement(elem));
}

var stats = mesh.GetColoringStatistics<Element, Node>();
Console.WriteLine($"Colors: {stats.NumberOfColors}, Avg group: {stats.AvgGroupSize}");
```

### ColoringStatistics

```csharp
public readonly struct ColoringStatistics
{
    public int ElementCount { get; }
    public int NumberOfColors { get; }
    public int MinGroupSize { get; }
    public int MaxGroupSize { get; }
    public double AvgGroupSize { get; }
}
```

---

## 15. Dual Graph Construction

```csharp
public DualGraph BuildDualGraph<TElement, TNode>(int minSharedNodes = 1);
public DualGraph BuildFaceNeighborGraph<TElement, TNode>();
public DualGraph BuildEdgeNeighborGraph<TElement, TNode>();
public DualGraph BuildVertexNeighborGraph<TElement, TNode>();
```

```csharp
// General dual graph (elements sharing ≥1 node are connected)
var dual = mesh.BuildDualGraph<Element, Node>(minSharedNodes: 1);

// Face-neighbor graph (for domain decomposition)
var faceDual = mesh.BuildFaceNeighborGraph<Element, Node>();

// Edge-neighbor or vertex-neighbor
var edgeDual = mesh.BuildEdgeNeighborGraph<Element, Node>();
var vertexDual = mesh.BuildVertexNeighborGraph<Element, Node>();
```

### DualGraph Class

```csharp
public sealed class DualGraph
{
    public List<List<int>> Adjacency { get; }
    public Dictionary<(int, int), int> SharedNodeCounts { get; }
    public int ElementCount { get; }
    public int EdgeCount { get; }

    public IReadOnlyList<int> GetNeighbors(int elementIndex);
    public int GetSharedNodeCount(int elem1, int elem2);
    public List<int> BreadthFirstSearch(int startElement);
    public List<List<int>> FindConnectedComponents();
    public int[] ComputeDistances(int sourceElement);
    public int ComputeDiameter();
}
```

```csharp
var dual = mesh.BuildDualGraph<Element, Node>();

// BFS on dual graph
var visited = dual.BreadthFirstSearch(0);

// Connected components
var components = dual.FindConnectedComponents();

// Graph diameter
int diameter = dual.ComputeDiameter();
```

---

## 16. Bandwidth Reduction

```csharp
public int[] ComputeCuthillMcKeeOrdering<TElement, TNode>(bool reverse = true);
public int ComputeBandwidth<TElement, TNode>();
public long ComputeProfile<TElement, TNode>();
public void ApplyNodePermutation<TElement, TNode>(int[] permutation);
```

```csharp
// Compute Reverse Cuthill-McKee ordering
int[] permutation = mesh.ComputeCuthillMcKeeOrdering<Element, Node>(reverse: true);

int bandwidthBefore = mesh.ComputeBandwidth<Element, Node>();
mesh.ApplyNodePermutation<Element, Node>(permutation);
int bandwidthAfter = mesh.ComputeBandwidth<Element, Node>();

Console.WriteLine($"Bandwidth reduced from {bandwidthBefore} to {bandwidthAfter}");
```

---

## 17. Sparse Matrix Patterns

```csharp
public (int[] RowPtr, int[] ColIndices) GetSparsityPatternCSR<TElement, TNode>(int dofsPerNode = 1);
public int GetNonZeroCount<TElement, TNode>(int dofsPerNode = 1);
public List<List<int>> GetCliques<TElement, TNode>();
public (int[] RowPtr, int[] ColumnIndices) ToCsr<TElement, TNode>();
public static Topology<TTypes> FromCsr<TElement, TNode>(int[] rowPtr, int[] columnIndices);
```

```csharp
// CSR sparsity pattern for FEA assembly
var (rowPtr, colIdx) = mesh.GetSparsityPatternCSR<Element, Node>(dofsPerNode: 3);
int nnz = mesh.GetNonZeroCount<Element, Node>(dofsPerNode: 3);

// Cliques (element node sets) for assembly
var cliques = mesh.GetCliques<Element, Node>();

// Export/import CSR format
var (rp, ci) = mesh.ToCsr<Element, Node>();
var imported = Topology<TypeMap<Node, Element>>.FromCsr<Element, Node>(rp, ci);
```

---

# Part IV: Smart Handles and Traversal

## 18. Smart Entity Handles

### 18.1 Creating Smart Handles

```csharp
public SmartEntity<TEntity> GetEntity<TEntity>(int index);
public List<SmartEntity<TEntity>> GetEntities<TEntity>(List<int> indices);
public List<SmartEntity<TEntity>> GetActiveEntities<TEntity>();
```

### 18.2 SmartEntity\<TEntity\> Record Struct

A fluent, object-oriented handle for navigating topology.

```csharp
public readonly record struct SmartEntity<TEntity> : IComparable<SmartEntity<TEntity>>
{
    public Topology<TTypes> Topology { get; init; }
    public int Index { get; init; }
    public bool IsValid { get; }
    public bool IsMarked { get; }
    public int Count { get; }

    // Data access
    public TData Data<TData>();
    public void SetData<TData>(TData value);

    // Navigation
    public List<SmartEntity<TRelated>> IncidentTo<TRelated>();
    public List<SmartEntity<TRelated>> Contains<TRelated>();
    public List<SmartEntity<TEntity>> Neighbors<TRelated>(bool sorted = true);
    public List<SmartEntity<TEntity>> DirectNeighbors<TRelated>(bool includeSelf = false, bool sorted = true);
    public List<(SmartEntity<TEntity> Entity, int SharedCount)> WeightedNeighbors<TRelated>();

    // Graph algorithms
    public Dictionary<SmartEntity<TEntity>, int> KHopNeighborhood<TRelated>(int k, int minShared = 1, bool includeSelf = true);
    public List<SmartEntity<TEntity>> EntitiesAtDistance<TRelated>(int k, int minShared = 1);
    public List<SmartEntity<TEntity>> BreadthFirstSearch(Action<SmartEntity<TEntity>, int>? visitor = null);
    public Dictionary<SmartEntity<TEntity>, int> BreadthFirstDistances();
    public Dictionary<SmartEntity<TEntity>, (double Distance, SmartEntity<TEntity> Predecessor)>
        DijkstraShortestPaths(Func<SmartEntity<TEntity>, SmartEntity<TEntity>, int, double> edgeWeight);

    // Modification
    public void MarkForRemoval();

    // Conversion
    public static implicit operator int(SmartEntity<TEntity> entity);
}
```

### 18.3 Usage

```csharp
// Get smart handles
var elem = mesh.GetEntity<Element>(0);
var allNodes = mesh.GetActiveEntities<Node>();

// Fluent navigation
Point coord = elem.Data<Point>();
var nodes = elem.Contains<Node>();
var neighbors = elem.Neighbors<Node>();

// Chain operations
var farNeighbors = elem
    .Neighbors<Node>(sorted: true)
    .SelectMany(n => n.Neighbors<Node>())
    .Distinct()
    .ToList();

// BFS from a handle
var visited = elem.BreadthFirstSearch((entity, depth) =>
{
    Console.WriteLine($"Element {entity.Index} at depth {depth}");
});

// Dijkstra
var paths = elem.DijkstraShortestPaths((from, to, edgeIdx) => 1.0);

// Implicit conversion to int
int idx = elem;  // Works seamlessly
```

---

## 19. Traversal and Components

```csharp
public IReadOnlyList<int> Traverse<TElement, TNode>(int startNode);
public IReadOnlyList<int> TraverseBreadthFirst<TElement, TNode>(int startNode);
public IReadOnlyList<IReadOnlyList<int>> FindComponents<TElement, TNode>();
```

```csharp
// Depth-first traversal from a node
var dfsOrder = mesh.Traverse<Element, Node>(startNodeIdx);

// Breadth-first traversal
var bfsOrder = mesh.TraverseBreadthFirst<Element, Node>(startNodeIdx);

// Find connected components
var components = mesh.FindComponents<Element, Node>();
Console.WriteLine($"Mesh has {components.Count} connected components");
```

---

# Part V: Operations and Performance

## 20. Batch Operations

```csharp
public void WithBatch(Action action);
public TResult WithBatch<TResult>(Func<TResult> func);
```

Batch operations acquire a single write lock for the entire block, providing 2–5x speedup for bulk operations.

```csharp
// Without batch: 10,000 lock/unlock cycles
for (int i = 0; i < 10_000; i++)
    mesh.Add<Node, Point>(points[i]);

// With batch: single lock acquisition (2-5x faster)
mesh.WithBatch(() =>
{
    for (int i = 0; i < 10_000; i++)
        mesh.Add<Node, Point>(points[i]);
});

// With return value
int count = mesh.WithBatch(() =>
{
    int added = 0;
    foreach (var pt in points)
    {
        mesh.Add<Node, Point>(pt);
        added++;
    }
    return added;
});
```

---

## 21. Compression and Optimization

### 21.1 Compress

```csharp
public void Compress(bool removeDuplicates = false, bool shrinkMemory = false, bool validate = false);
public async Task CompressAsync(bool removeDuplicates = false, bool shrinkMemory = false,
    bool validate = false, CancellationToken cancellationToken = default);
public void Clear();
```

```csharp
// Simple: apply pending removals and renumber
mesh.Compress();

// Full optimization
mesh.Compress(removeDuplicates: true, shrinkMemory: true, validate: true);

// Async version
await mesh.CompressAsync(removeDuplicates: true, shrinkMemory: true);

// Clear everything
mesh.Clear();
```

### 21.2 Memory Management

```csharp
public void ConfigureType<TEntity>(int parallelizationThreshold, int? reserveCapacity = null);
public void Reserve<TElement, TNode>(int capacity);
public void ShrinkToFit();
public long EstimateMemoryUsage();
```

```csharp
// Pre-allocate capacity
mesh.Reserve<Element, Node>(expectedCount: 100_000);
mesh.ConfigureType<Element>(parallelizationThreshold: 5000, reserveCapacity: 100_000);

// Reclaim excess memory
mesh.ShrinkToFit();

// Estimate memory usage
long bytes = mesh.EstimateMemoryUsage();
Console.WriteLine($"Memory: {bytes / 1024.0 / 1024.0:F1} MB");
```

### 21.3 Statistics

```csharp
public TopologyStats GetStatistics();
public ElementStatistics GetElementStatistics<TElement, TNode>();
```

### TopologyStats

```csharp
public sealed class TopologyStats
{
    public IReadOnlyDictionary<Type, int> EntityCounts { get; }
    public IReadOnlyDictionary<Type, int> DataCounts { get; }
    public IReadOnlyList<Type> TypesWithSymmetry { get; }
    public int TotalEntities { get; }
}
```

### ElementStatistics

```csharp
public readonly struct ElementStatistics
{
    public int ElementCount { get; }
    public int MinNodesPerElement { get; }
    public int MaxNodesPerElement { get; }
    public double AvgNodesPerElement { get; }
    public IReadOnlyDictionary<int, int> NodesPerElementDistribution { get; }
}
```

---

## 22. Mesh Merging and Extraction

```csharp
public int Merge<TElement, TNode>(Topology<TTypes> other);

public (Topology<TTypes> Subtopology, int[] NodeMap, int[] ElementMap)
    ExtractSubstructure<TElement, TNode>(IEnumerable<int> elementIndices);

public (Topology<TTypes> NewTopology, int[] ElementMapping, int[] NodeMapping)
    CloneWhere<TElement, TNode>(Func<int, bool> predicate, bool includeOrphanNodes = false);

public (Topology<TTypes> NewTopology, int[] ElementMapping, int[] NodeMapping)
    ExtractRegion<TElement, TNode>(IEnumerable<int> elementIndices);

public (Topology<TTypes> NewTopology, int[] ElementMapping, int[] NodeMapping)
    ExtractByBoundingBox<TElement, TNode>(
        (double X, double Y, double Z) minBound,
        (double X, double Y, double Z) maxBound,
        Func<int, (double X, double Y, double Z)> getCoord,
        bool allNodesInside = false);
```

```csharp
// Merge another topology
int nodeOffset = mesh.Merge<Element, Node>(otherMesh);

// Extract sub-mesh by element indices
var (subMesh, nodeMap, elemMap) = mesh.ExtractSubstructure<Element, Node>([0, 1, 5, 10]);

// Clone with predicate
var (filtered, eMap, nMap) = mesh.CloneWhere<Element, Node>(
    elemIdx => mesh.Get<Element, Material>(elemIdx).E > 100e9);

// Extract by bounding box
var (region, eMapping, nMapping) = mesh.ExtractByBoundingBox<Element, Node>(
    minBound: (0, 0, 0),
    maxBound: (1, 1, 1),
    getCoord: nodeIdx =>
    {
        var p = mesh.Get<Node, Point>(nodeIdx);
        return (p.X, p.Y, p.Z);
    });
```

---

## 23. Serialization

```csharp
public string ToJson(JsonSerializerOptions? options = null);
public void SaveToFile(string path, JsonSerializerOptions? options = null);
public static Topology<TTypes> FromJson(string json, JsonSerializerOptions? options = null);
public static Topology<TTypes> LoadFromFile(string path, JsonSerializerOptions? options = null);
```

```csharp
// Save to JSON string
string json = mesh.ToJson();

// Save to file
mesh.SaveToFile("mesh.json");

// Load from file
var loaded = Topology<TypeMap<Node, Element>>.LoadFromFile("mesh.json");
```

---

## 24. Performance Configuration

### ParallelConfig Static Class

```csharp
public static class ParallelConfig
{
    // Thread control
    public static int MaxDegreeOfParallelism { get; set; }
    public static ParallelOptions Options { get; }
    public static int ProcessorCount { get; }

    // GPU
    public static bool EnableGPU { get; set; }
    public static bool UseGPU { get; }
    public static bool IsGPUAvailable { get; }
    public static bool SkipGPUCheck { get; set; }

    // MKL
    public static int MKLNumThreads { get; set; }
    public static int? MKLCurrentThreads { get; }
    public static bool IsMKLAvailable { get; }

    // Debug
    public static bool EnableDebugOutput { get; set; }

    // Convenience
    public static void SetAllThreads(int numThreads);
    public static void Reset();
    public static void Cleanup();
    public static string GetSummary();
}
```

```csharp
ParallelConfig.MaxDegreeOfParallelism = 8;
ParallelConfig.SetAllThreads(4);
Console.WriteLine(ParallelConfig.GetSummary());  // "CPU=4/8, GPU=true, MKL=4"
ParallelConfig.Reset();
```

---

# Part VI: Supporting Types

## 25. O2M — One-to-Many Sparse Structure

`O2M` is a sparse adjacency structure representing one-to-many relationships. It is the algebraic foundation for connectivity operations.

```csharp
public sealed class O2M : IComparable<O2M>, IEquatable<O2M>, ICloneable
```

### 25.1 Construction

```csharp
public O2M();
public O2M(int reservedCapacity);
public O2M(List<List<int>> adjacenciesList);

// Factory methods
public static O2M FromCsr(int[] rowPointers, int[] columnIndices);
public static O2M FromBooleanMatrix(bool[,] matrix);
public static O2M GetRandomO2M(int elementCount, int nodeCount, double density, int? seed = null);

// Implicit conversions
public static implicit operator O2M(List<List<int>> nodes);
public static implicit operator O2M(List<int> elements);
```

### 25.2 Access and Query

```csharp
public int Count { get; }
public ReadOnlySpan<int> this[int rowIndex] { get; }
public int this[int rowIndex, int columnIndex] { get; }

public int GetMaxNode();
public long GetTotalEdgeCount();
public (int MinDegree, int MaxDegree, double AvgDegree, long TotalEdges) GetStatistics();
public bool IsValid();
public bool IsSorted();
public string? ValidateStrict();
public bool IsAcyclic();
public List<int> GetTopOrder();
public List<int> GetSortOrder();
public List<int> GetDuplicates();
```

### 25.3 Modification

```csharp
public int AppendElement(List<int> nodes);
public int AppendElementCopy(ReadOnlySpan<int> nodes);
public void AppendElements(params List<int>[] nodes);
public void AppendNodeToElement(int elementIndex, int nodeValue);
public bool RemoveNodeFromElement(int elementIndex, int nodeValue);
public void ClearElement(int elementIndex);
public void ReplaceElement(int elementIndex, List<int> newNodes);
public void ReplaceElementCopy(int elementIndex, ReadOnlySpan<int> newNodes);
public void ClearAll();
public void Reserve(int reservedCapacity);
public void ShrinkToFit();
```

### 25.4 Operators (Algebraic)

```csharp
// Matrix multiplication (element-wise AND)
public static O2M operator *(O2M left, O2M right);

// Set operations
public static O2M operator +(O2M left, O2M right);   // Union
public static O2M operator &(O2M left, O2M right);   // Intersection
public static O2M operator -(O2M left, O2M right);   // Difference
public static O2M operator ^(O2M left, O2M right);   // Symmetric difference
public static O2M operator |(O2M left, O2M right);   // Element-wise OR
```

```csharp
O2M A = elemToNode;
O2M At = A.Transpose();
O2M AAt = A * At;     // Element-to-element via shared nodes

O2M union = A + B;
O2M intersect = A & B;
O2M diff = A - B;
```

### 25.5 Transpose

```csharp
public O2M Transpose();
public O2M Transpose(int maxNodeCap);
public O2M TransposeStrict();
```

### 25.6 Permutation and Compression

```csharp
public void PermuteElements(List<int> oldToNewElementMap);
public unsafe void PermuteNodes(List<int> oldToNewNodeMap);
public void CompressElements(List<int> newToOldElementMap);
public void RearrangeAfterRenumbering(List<int> newToOldElementMap, List<int> oldToNewNodeMap);
```

### 25.7 Conversion

```csharp
public bool[,] ToBooleanMatrix();
public List<List<int>> ToAdjacencyLists();
public string ToEpsString();
public object Clone();
```

### 25.8 Comparison

```csharp
public int CompareTo(O2M? other);
public bool Equals(O2M? other);
public bool IsPermutationOf(O2M? other);
public int FullContentHashCode();
```

### 25.9 Graph Algorithms on O2M

```csharp
public List<int> BreadthFirstSearch(int startElement, Action<int, int>? visitor = null);
public List<int> BreadthFirstSearch(int startElement, O2M? transpose, Action<int, int>? visitor = null);
public Dictionary<int, int> BreadthFirstDistances(int startElement);
public Dictionary<int, int> BreadthFirstDistances(int startElement, O2M? transpose);
public Dictionary<int, (double Distance, int Predecessor)> DijkstraShortestPaths(
    int startElement, Func<int, int, int, double> edgeWeight);
public Dictionary<int, (double Distance, int Predecessor)> DijkstraShortestPaths(
    int startElement, O2M? transpose, Func<int, int, int, double> edgeWeight);
public List<List<int>> GetCliques();
```

### 25.10 Static Analysis Methods

```csharp
public static List<List<int>> GetNodePositions(O2M nodesFromElement, O2M elementsFromNode);
public static List<List<int>> GetElementPositions(O2M nodesFromElement, O2M elementsFromNode);
public static List<List<int>> GetCliques(O2M nodesFromElement, O2M elementsFromNode);
public static List<List<int>> GetCliquesStrict(O2M nodesFromElement, O2M elementsFromNode);
```

---

## 26. M2M — Thread-Safe Many-to-Many

`M2M` wraps `O2M` with `ReaderWriterLockSlim` for thread-safe access to element–node relationships.

```csharp
public sealed class M2M : IComparable<M2M>, IEquatable<M2M>, IDisposable
```

### 26.1 Construction

```csharp
public M2M();
public M2M(int reservedCapacity);
public M2M(List<List<int>> adjacencies);
public M2M(O2M nodesFromElement);

public static M2M FromCsr(int[] rowPointers, int[] columnIndices);
public static M2M FromBooleanMatrix(bool[,] matrix);
public static M2M CreateRandom(int elementCount, int nodeCount, double density, int? seed = null);
```

### 26.2 Element and Node Access

```csharp
public IReadOnlyList<int> GetNodesForElement(int elementIndex);
public void WithNodesSpan(int elementIndex, Action<ReadOnlySpan<int>> action);
public TResult WithNodesSpan<TResult>(int elementIndex, Func<ReadOnlySpan<int>, TResult> func);

public int GetTransposeNodeCount();
public int GetElementCountForNode(int node);
public List<int> GetElementsForNode(int node);
public void WithElementsForNode(int node, Action<IReadOnlyList<int>> action);
public void WithElementsForNodeSpan(int node, ReadOnlySpanAction<int> action);
```

### 26.3 Queries

```csharp
public List<int> GetElementsWithNodes(List<int> nodes);
public List<int> GetElementsContainingAnyNode(List<int> nodes);
public List<int> GetElementsFromNodes(List<int> nodes);
public List<int> GetElementNeighbors(int element, bool sorted = true);
public List<int> GetNodeNeighbors(int node, bool sorted = true);
public bool HasElement(int elementIndex);
public bool HasNode(int nodeIndex);
public bool ElementContainsNode(int elementIndex, int nodeIndex);
public int GetMaxNode();
public bool IsValid();
```

### 26.4 Graph Operations

```csharp
public O2M GetElementsToElements();
public M2M GetElementToElementGraph();
public O2M GetNodesToNodes();
public M2M GetNodeToNodeGraph();
public List<int> GetTopologicalOrder();
public bool IsAcyclic();
public List<int> GetSortOrder();
public List<int> GetDuplicates();
public bool IsPermutationOf(M2M? other);
```

### 26.5 Modification

```csharp
public int AppendElement(List<int> nodes);
public void AppendElements(params List<int>[] nodes);
public void AppendNodeToElement(int elementIndex, int nodeValue);
public bool RemoveNodeFromElement(int elementIndex, int nodeValue);
public void ClearElement(int elementIndex);
public void ReplaceElement(int elementIndex, List<int> newNodes);
public void ClearAll();
public void Reserve(int capacity);
public void ShrinkToFit();
```

### 26.6 Permutation and Batch

```csharp
public void CompressElements(List<int> newToOldElementMap);
public void PermuteElements(List<int> oldToNewElementMap);
public void PermuteNodes(List<int> oldToNewNodeMap);
public void RearrangeAfterRenumbering(List<int> newToOldElementMap, List<int> oldToNewNodeMap);
public IDisposable BeginBatchUpdate();
public void Synchronize();
public void EnsurePositionCaches();
```

### 26.7 Graph Algorithms

```csharp
public List<int> BreadthFirstSearch(int startElement, Action<int, int>? visitor = null);
public Dictionary<int, int> BreadthFirstDistances(int startElement);
public Dictionary<int, (double Distance, int Predecessor)> DijkstraShortestPaths(
    int startElement, Func<int, int, int, double> edgeWeight);
```

### 26.8 Conversion and Comparison

```csharp
public bool[,] ToBooleanMatrix();
public object Clone();
public M2M CloneTyped();
public int CompareTo(M2M? other);
public bool Equals(M2M? other);
public void Dispose();
```

---

## 27. MM2M — Multi-Type Manager

`MM2M` manages an N×N matrix of `M2M` blocks for multi-type topologies.

```csharp
public sealed class MM2M : IDisposable
```

### 27.1 Construction and Access

```csharp
public MM2M(int numberOfTypes);

public int NumberOfTypes { get; }
public int Count { get; }
public M2M this[int elementType, int nodeType] { get; set; }
```

### 27.2 Block Operations

```csharp
public void WithBlock(int elementType, int nodeType, Action<M2M> action);
public TResult WithBlock<TResult>(int elementType, int nodeType, Func<M2M, TResult> func);
```

### 27.3 Cross-Type Queries

```csharp
public List<(int ElemType, int Elem)> GetAllElements(int nodeType, int node, bool sorted = true);
public int GetNumberOfNodes(int elementType, int element, int nodeType);
public bool TryGetNumberOfNodes(int elementType, int element, int nodeType, out int count);
public int GetNumberOfElements(int nodeType, int node, int elementType);
```

### 27.4 Multi-Type Graph Algorithms

```csharp
public List<(int ElementType, int Element)> BreadthFirstSearchMultiType(
    int startElementType, int startElement, Action<int, int, int>? visitor = null);

public Dictionary<(int ElementType, int Element), int> BreadthFirstDistancesMultiType(
    int startElementType, int startElement);

public Dictionary<(int ElementType, int Element), (double Distance, (int PredType, int PredElem))>
    DijkstraShortestPathsMultiType(int startElementType, int startElement,
    Func<int, int, int, int, int, int, double> edgeWeight);

public static List<(int ElementType, int Element)>? ReconstructPathMultiType(
    Dictionary<(int ElementType, int Element), (double Distance, (int PredType, int PredElem))> dijkstraResult,
    int targetElementType, int targetElement);
```

---

## 28. Utility Functions

### Utils Static Class

```csharp
public static class Utils
{
    // List operations
    public static void SortUnique<T>(this List<T> list) where T : IComparable<T>;
    public static void InsertSorted(this List<int> sortedList, int value);
    public static bool RemoveSorted(this List<int> sortedList, int value);
    public static bool ContainsSorted(List<int> sortedList, int value);
    public static int BinarySearch(List<int> sortedList, int value);

    // Sorted set operations
    public static List<int> UnionSorted(List<int> a, List<int> b);
    public static List<int> IntersectSorted(List<int> a, List<int> b);
    public static List<int> DifferenceSorted(List<int> a, List<int> b);

    // Comparison
    public static int Compare<T>(List<T> first, List<T> second) where T : IComparable<T>;
    public static bool AreEqual<T>(List<T> first, List<T> second) where T : IComparable<T>;

    // Aggregation
    public static int Min(List<int> list);
    public static int Max(List<int> list);
    public static long Sum(List<int> list);
    public static int IndexOfMin(List<int> list);
    public static int IndexOfMax(List<int> list);

    // Copy
    public static List<T> Copy<T>(List<T> source);
    public static List<List<T>> DeepCopy<T>(List<List<T>> source);
    public static List<T> ToList<T>(HashSet<T> set);
    public static List<T> ToSortedList<T>(HashSet<T> set) where T : IComparable<T>;

    // Helpers
    public static List<int> Range(int count);
    public static bool AreIndicesValid(List<int> indices, int collectionSize);
    public static bool AreAllNonNegative(List<int> list);
    public static List<T> GetItemsAtIndices<T>(this List<T> list, List<int> indices);

    // Node renumbering
    public static (List<int> newNodesFromOld, List<int> oldNodesFromNew)
        GetNodeMapsFromKillList(List<int> killList);
}
```

### ReadOnlySpanAction Delegate

```csharp
public delegate void ReadOnlySpanAction<T>(ReadOnlySpan<T> span);
```

---

# Appendices

## Appendix A: API Quick Reference

### Entity Operations

| Method | Description | Complexity |
|--------|-------------|------------|
| `Add<TNode>()` | Add standalone entity | O(1) amortized |
| `Add<TElement, TNode>(...)` | Add connected entity | O(1) amortized |
| `AddUnique<TElement, TNode>(...)` | Add with deduplication | O(m + log α) |
| `AddRange<TNode, TData>(...)` | Bulk add with data | O(n) |
| `AddRangeParallel<TElement, TNode>(...)` | Parallel bulk add | O(n/p) |
| `Remove<TEntity>(idx)` | Mark for removal | O(1) |
| `Compress()` | Apply removals | O(n) |

### Connectivity Queries

| Method | Description | Complexity |
|--------|-------------|------------|
| `NodesOf<TElement, TNode>(elem)` | Element → nodes | O(1) |
| `ElementsAt<TElement, TNode>(node)` | Node → elements | O(1) cached |
| `Neighbors<TElement, TNode>(elem)` | Neighbor elements | O(M·K) |
| `CountRelated<TEntity, TRelated>(idx)` | Count connections | O(1) |
| `CountIncident<TElement, TNode>(node)` | Count incident elements | O(1) |

### Data Operations

| Method | Description | Complexity |
|--------|-------------|------------|
| `Get<TEntity, TData>(idx)` | Get data | O(1) |
| `Set<TEntity, TData>(idx, val)` | Set data | O(1) |
| `All<TEntity, TData>()` | Get all data | O(n) |
| `ForEach<TEntity, TData>(action)` | Iterate with data | O(n) |
| `ParallelForEach<TEntity, TData>(action)` | Parallel iterate | O(n/p) |

### Advanced Operations

| Method | Description |
|--------|-------------|
| `FindConnectedComponents<>()` | Connected components |
| `ComputeElementColoring<>()` | Graph coloring |
| `ComputeCuthillMcKeeOrdering<>()` | Bandwidth reduction |
| `BuildDualGraph<>()` | Dual graph |
| `DiscoverSubEntities<>()` | Sub-entity extraction |
| `FindBoundaryNodes<>()` | Boundary detection |
| `GetSparsityPatternCSR<>()` | Sparse matrix pattern |
| `BreadthFirstSearch<>()` | BFS traversal |
| `DijkstraShortestPaths<>()` | Shortest paths |

---

## Appendix B: Performance Characteristics

### Throughput Benchmarks (Typical Hardware)

| Operation | Rate | Notes |
|-----------|------|-------|
| Serial Add | ~120K elem/sec | Individual additions |
| Batch Add | ~430K elem/sec | 3.6x faster than serial |
| Parallel Add | ~690K elem/sec | 5.8x faster (8 cores) |
| NodesOf query | ~12.5M/sec | 80 ns per query |
| ElementsAt query | ~8.3M/sec | Cached transpose |
| Neighbors query | ~320–420K/sec | Depends on connectivity |

### Memory Usage (Typical 100K Tets, 25K Nodes)

| Component | Size |
|-----------|------|
| Nodes | ~600 KB |
| Elements | ~4 MB |
| Transpose cache | ~4 MB |
| Canonical index | ~4 MB |
| **Total** | **~13 MB** |

### Parallel Speedup (8 Cores)

| Operation | Speedup |
|-----------|---------|
| AddRangeParallel | 5.8x |
| Transpose | 6.2x |
| ParallelForEach | 7.1x |

---

## Appendix C: Common Patterns

### FEA Assembly Loop

```csharp
mesh.ForEach<Element, Material>((elemIdx, material) =>
{
    mesh.WithNodesSpan<Element, Node>(elemIdx, nodes =>
    {
        Span<Point> coords = stackalloc Point[nodes.Length];
        for (int i = 0; i < nodes.Length; i++)
            coords[i] = mesh.Get<Node, Point>(nodes[i]);

        var Ke = ComputeElementStiffness(coords, material);
        AssembleGlobal(Ke, nodes);
    });
});
```

### Edge Extraction with Deduplication

```csharp
using var mesh = Topology.New<Node, Edge, Element>();
mesh.WithSymmetry<Edge>(Symmetry.Dihedral(2));

// Add nodes and elements...

var triEdges = SubEntityDefinition.FromEdges((0,1), (1,2), (2,0));
mesh.DiscoverSubEntities<Element, Edge, Node>(triEdges);
```

### Race-Free Parallel Assembly

```csharp
var groups = mesh.GetColorGroups<Element, Node>();

foreach (var group in groups)
{
    Parallel.ForEach(group, elem =>
    {
        mesh.WithNodesSpan<Element, Node>(elem, nodes =>
        {
            var Ke = ComputeStiffness(nodes);
            AssembleToGlobal(Ke, nodes);  // No race conditions within a color group
        });
    });
}
```

### Bandwidth Reduction before Assembly

```csharp
int[] perm = mesh.ComputeCuthillMcKeeOrdering<Element, Node>(reverse: true);
mesh.ApplyNodePermutation<Element, Node>(perm);

var (rowPtr, colIdx) = mesh.GetSparsityPatternCSR<Element, Node>(dofsPerNode: 3);
// Use rowPtr, colIdx for sparse solver
```

### Mesh Refinement

```csharp
mesh.WithBatch(() =>
{
    foreach (int elem in elementsToRefine)
    {
        var nodes = mesh.NodesOf<Element, Node>(elem);
        // Create midpoint nodes
        int mid01 = mesh.Add<Node, Point>(Midpoint(nodes[0], nodes[1]));
        int mid12 = mesh.Add<Node, Point>(Midpoint(nodes[1], nodes[2]));
        int mid20 = mesh.Add<Node, Point>(Midpoint(nodes[2], nodes[0]));

        // Replace original with center triangle
        mesh.ReplaceElementNodes<Element, Node>(elem, mid01, mid12, mid20);

        // Add corner triangles
        mesh.Add<Element, Node>(nodes[0], mid01, mid20);
        mesh.Add<Element, Node>(mid01, nodes[1], mid12);
        mesh.Add<Element, Node>(mid20, mid12, nodes[2]);
    }
});
```

### Domain Decomposition

```csharp
var dual = mesh.BuildFaceNeighborGraph<Element, Node>();
var components = dual.FindConnectedComponents();

// Or extract sub-meshes by region
var (subMesh, nodeMap, elemMap) = mesh.ExtractByBoundingBox<Element, Node>(
    minBound: (0, 0, 0), maxBound: (0.5, 0.5, 0.5),
    getCoord: n => { var p = mesh.Get<Node, Point>(n); return (p.X, p.Y, p.Z); });
```

---

*This document covers the complete public API of the Topology library. All types, methods, and properties listed are part of the public contract and available for use in applications.*
