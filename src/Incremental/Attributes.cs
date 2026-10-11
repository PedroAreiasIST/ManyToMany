// Attributes.cs — per-element and per-node data versioned with the mesh.
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
//
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
