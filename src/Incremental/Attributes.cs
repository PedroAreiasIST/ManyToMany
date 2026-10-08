// Attributes.cs — per-element and per-node data versioned with the mesh.
//
// An attribute is a persistent array registered with an IncrementalMesh: it grows
// with the element or node count, is captured by Commit and restored by Rollback,
// so data such as refinement tags, coordinates or history variables cannot drift
// from the topology they describe.

namespace Step2;

public enum AttributeKind
{
    Element,
    Node
}

internal interface IAttribute
{
    AttributeKind Kind { get; }
    void Resize(int length);
    object Freeze();
    void Restore(object view);
    void PrepareWrite();
}

public sealed class Attribute<T> : IAttribute
{
    internal readonly PArray<T> Data;
    public readonly int Index;
    public AttributeKind Kind { get; }

    internal Attribute(PArray<T> data, AttributeKind kind, int index)
    {
        Data = data;
        Kind = kind;
        Index = index;
    }

    public T this[int i] => Data[i];
    public void Set(int i, T value) => Data.Set(i, value);
    public int Length => Data.Length;

    void IAttribute.Resize(int length) => Data.Resize(length);
    object IAttribute.Freeze() => Data.Freeze();
    void IAttribute.Restore(object view) => Data.Restore((PView<T>)view);
    void IAttribute.PrepareWrite() => Data.PrepareWrite();
}
