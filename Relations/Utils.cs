using System.Buffers;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Diagnostics;
using System.Reflection;

namespace Numerical;

#region TypeMap - Compile-time type-to-integer mapping

/// <summary>
///     Compile-time type-to-integer mapping for use with Topology.
/// </summary>
public interface ITypeMap
{
    int Count { get; }

    /// <summary>
    ///     Gets the index of type T in this map.
    /// </summary>
    /// <exception cref="KeyNotFoundException">Thrown if T is not in this TypeMap.</exception>
    int IndexOf<T>()
    {
        if (TryIndexOf<T>(out var index)) return index;
        throw new KeyNotFoundException($"Type {typeof(T).Name} is not in this TypeMap.");
    }

    /// <summary>
    ///     Tries to get the index of type T in this map.
    /// </summary>
    /// <param name="index">The index if found; -1 otherwise.</param>
    /// <returns>True if T is in this TypeMap; false otherwise.</returns>
    bool TryIndexOf<T>(out int index);
}

/// <summary>Type map for 2 types.</summary>
public sealed class TypeMap<T0, T1> : ITypeMap
{
    public int Count => 2;


    public bool TryIndexOf<T>(out int index)
    {
        if (typeof(T) == typeof(T0))
        {
            index = 0;
            return true;
        }

        if (typeof(T) == typeof(T1))
        {
            index = 1;
            return true;
        }

        index = -1;
        return false;
    }
}

/// <summary>Type map for 3 types.</summary>
public sealed class TypeMap<T0, T1, T2> : ITypeMap
{
    public int Count => 3;


    public bool TryIndexOf<T>(out int index)
    {
        if (typeof(T) == typeof(T0))
        {
            index = 0;
            return true;
        }

        if (typeof(T) == typeof(T1))
        {
            index = 1;
            return true;
        }

        if (typeof(T) == typeof(T2))
        {
            index = 2;
            return true;
        }

        index = -1;
        return false;
    }
}

/// <summary>Type map for 4 types.</summary>
public sealed class TypeMap<T0, T1, T2, T3> : ITypeMap
{
    public int Count => 4;


    public bool TryIndexOf<T>(out int index)
    {
        if (typeof(T) == typeof(T0))
        {
            index = 0;
            return true;
        }

        if (typeof(T) == typeof(T1))
        {
            index = 1;
            return true;
        }

        if (typeof(T) == typeof(T2))
        {
            index = 2;
            return true;
        }

        if (typeof(T) == typeof(T3))
        {
            index = 3;
            return true;
        }

        index = -1;
        return false;
    }
}

/// <summary>Type map for 5 types.</summary>
public sealed class TypeMap<T0, T1, T2, T3, T4> : ITypeMap
{
    public int Count => 5;


    public bool TryIndexOf<T>(out int index)
    {
        if (typeof(T) == typeof(T0))
        {
            index = 0;
            return true;
        }

        if (typeof(T) == typeof(T1))
        {
            index = 1;
            return true;
        }

        if (typeof(T) == typeof(T2))
        {
            index = 2;
            return true;
        }

        if (typeof(T) == typeof(T3))
        {
            index = 3;
            return true;
        }

        if (typeof(T) == typeof(T4))
        {
            index = 4;
            return true;
        }

        index = -1;
        return false;
    }
}

/// <summary>Type map for 6 types.</summary>
public sealed class TypeMap<T0, T1, T2, T3, T4, T5> : ITypeMap
{
    public int Count => 6;


    public bool TryIndexOf<T>(out int index)
    {
        if (typeof(T) == typeof(T0))
        {
            index = 0;
            return true;
        }

        if (typeof(T) == typeof(T1))
        {
            index = 1;
            return true;
        }

        if (typeof(T) == typeof(T2))
        {
            index = 2;
            return true;
        }

        if (typeof(T) == typeof(T3))
        {
            index = 3;
            return true;
        }

        if (typeof(T) == typeof(T4))
        {
            index = 4;
            return true;
        }

        if (typeof(T) == typeof(T5))
        {
            index = 5;
            return true;
        }

        index = -1;
        return false;
    }
}

/// <summary>Type map for 7 types.</summary>
public sealed class TypeMap<T0, T1, T2, T3, T4, T5, T6> : ITypeMap
{
    public int Count => 7;


    public bool TryIndexOf<T>(out int index)
    {
        if (typeof(T) == typeof(T0))
        {
            index = 0;
            return true;
        }

        if (typeof(T) == typeof(T1))
        {
            index = 1;
            return true;
        }

        if (typeof(T) == typeof(T2))
        {
            index = 2;
            return true;
        }

        if (typeof(T) == typeof(T3))
        {
            index = 3;
            return true;
        }

        if (typeof(T) == typeof(T4))
        {
            index = 4;
            return true;
        }

        if (typeof(T) == typeof(T5))
        {
            index = 5;
            return true;
        }

        if (typeof(T) == typeof(T6))
        {
            index = 6;
            return true;
        }

        index = -1;
        return false;
    }
}

/// <summary>Type map for 8 types.</summary>
public sealed class TypeMap<T0, T1, T2, T3, T4, T5, T6, T7> : ITypeMap
{
    public int Count => 8;


    public bool TryIndexOf<T>(out int index)
    {
        if (typeof(T) == typeof(T0))
        {
            index = 0;
            return true;
        }

        if (typeof(T) == typeof(T1))
        {
            index = 1;
            return true;
        }

        if (typeof(T) == typeof(T2))
        {
            index = 2;
            return true;
        }

        if (typeof(T) == typeof(T3))
        {
            index = 3;
            return true;
        }

        if (typeof(T) == typeof(T4))
        {
            index = 4;
            return true;
        }

        if (typeof(T) == typeof(T5))
        {
            index = 5;
            return true;
        }

        if (typeof(T) == typeof(T6))
        {
            index = 6;
            return true;
        }

        if (typeof(T) == typeof(T7))
        {
            index = 7;
            return true;
        }

        index = -1;
        return false;
    }
}

/// <summary>Type map for 9 types.</summary>
public sealed class TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8> : ITypeMap
{
    public int Count => 9;


    public bool TryIndexOf<T>(out int index)
    {
        if (typeof(T) == typeof(T0))
        {
            index = 0;
            return true;
        }

        if (typeof(T) == typeof(T1))
        {
            index = 1;
            return true;
        }

        if (typeof(T) == typeof(T2))
        {
            index = 2;
            return true;
        }

        if (typeof(T) == typeof(T3))
        {
            index = 3;
            return true;
        }

        if (typeof(T) == typeof(T4))
        {
            index = 4;
            return true;
        }

        if (typeof(T) == typeof(T5))
        {
            index = 5;
            return true;
        }

        if (typeof(T) == typeof(T6))
        {
            index = 6;
            return true;
        }

        if (typeof(T) == typeof(T7))
        {
            index = 7;
            return true;
        }

        if (typeof(T) == typeof(T8))
        {
            index = 8;
            return true;
        }

        index = -1;
        return false;
    }
}

/// <summary>Type map for 10 types.</summary>
public sealed class TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9> : ITypeMap
{
    public int Count => 10;


    public bool TryIndexOf<T>(out int index)
    {
        if (typeof(T) == typeof(T0))
        {
            index = 0;
            return true;
        }

        if (typeof(T) == typeof(T1))
        {
            index = 1;
            return true;
        }

        if (typeof(T) == typeof(T2))
        {
            index = 2;
            return true;
        }

        if (typeof(T) == typeof(T3))
        {
            index = 3;
            return true;
        }

        if (typeof(T) == typeof(T4))
        {
            index = 4;
            return true;
        }

        if (typeof(T) == typeof(T5))
        {
            index = 5;
            return true;
        }

        if (typeof(T) == typeof(T6))
        {
            index = 6;
            return true;
        }

        if (typeof(T) == typeof(T7))
        {
            index = 7;
            return true;
        }

        if (typeof(T) == typeof(T8))
        {
            index = 8;
            return true;
        }

        if (typeof(T) == typeof(T9))
        {
            index = 9;
            return true;
        }

        index = -1;
        return false;
    }
}

/// <summary>Type map for 11 types.</summary>
public sealed class TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10> : ITypeMap
{
    public int Count => 11;


    public bool TryIndexOf<T>(out int index)
    {
        if (typeof(T) == typeof(T0))
        {
            index = 0;
            return true;
        }

        if (typeof(T) == typeof(T1))
        {
            index = 1;
            return true;
        }

        if (typeof(T) == typeof(T2))
        {
            index = 2;
            return true;
        }

        if (typeof(T) == typeof(T3))
        {
            index = 3;
            return true;
        }

        if (typeof(T) == typeof(T4))
        {
            index = 4;
            return true;
        }

        if (typeof(T) == typeof(T5))
        {
            index = 5;
            return true;
        }

        if (typeof(T) == typeof(T6))
        {
            index = 6;
            return true;
        }

        if (typeof(T) == typeof(T7))
        {
            index = 7;
            return true;
        }

        if (typeof(T) == typeof(T8))
        {
            index = 8;
            return true;
        }

        if (typeof(T) == typeof(T9))
        {
            index = 9;
            return true;
        }

        if (typeof(T) == typeof(T10))
        {
            index = 10;
            return true;
        }

        index = -1;
        return false;
    }
}

/// <summary>Type map for 12 types.</summary>
public sealed class TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11> : ITypeMap
{
    public int Count => 12;


    public bool TryIndexOf<T>(out int index)
    {
        if (typeof(T) == typeof(T0))
        {
            index = 0;
            return true;
        }

        if (typeof(T) == typeof(T1))
        {
            index = 1;
            return true;
        }

        if (typeof(T) == typeof(T2))
        {
            index = 2;
            return true;
        }

        if (typeof(T) == typeof(T3))
        {
            index = 3;
            return true;
        }

        if (typeof(T) == typeof(T4))
        {
            index = 4;
            return true;
        }

        if (typeof(T) == typeof(T5))
        {
            index = 5;
            return true;
        }

        if (typeof(T) == typeof(T6))
        {
            index = 6;
            return true;
        }

        if (typeof(T) == typeof(T7))
        {
            index = 7;
            return true;
        }

        if (typeof(T) == typeof(T8))
        {
            index = 8;
            return true;
        }

        if (typeof(T) == typeof(T9))
        {
            index = 9;
            return true;
        }

        if (typeof(T) == typeof(T10))
        {
            index = 10;
            return true;
        }

        if (typeof(T) == typeof(T11))
        {
            index = 11;
            return true;
        }

        index = -1;
        return false;
    }
}

/// <summary>Type map for 13 types.</summary>
public sealed class TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12> : ITypeMap
{
    public int Count => 13;


    public bool TryIndexOf<T>(out int index)
    {
        if (typeof(T) == typeof(T0))
        {
            index = 0;
            return true;
        }

        if (typeof(T) == typeof(T1))
        {
            index = 1;
            return true;
        }

        if (typeof(T) == typeof(T2))
        {
            index = 2;
            return true;
        }

        if (typeof(T) == typeof(T3))
        {
            index = 3;
            return true;
        }

        if (typeof(T) == typeof(T4))
        {
            index = 4;
            return true;
        }

        if (typeof(T) == typeof(T5))
        {
            index = 5;
            return true;
        }

        if (typeof(T) == typeof(T6))
        {
            index = 6;
            return true;
        }

        if (typeof(T) == typeof(T7))
        {
            index = 7;
            return true;
        }

        if (typeof(T) == typeof(T8))
        {
            index = 8;
            return true;
        }

        if (typeof(T) == typeof(T9))
        {
            index = 9;
            return true;
        }

        if (typeof(T) == typeof(T10))
        {
            index = 10;
            return true;
        }

        if (typeof(T) == typeof(T11))
        {
            index = 11;
            return true;
        }

        if (typeof(T) == typeof(T12))
        {
            index = 12;
            return true;
        }

        index = -1;
        return false;
    }
}

/// <summary>Type map for 14 types.</summary>
public sealed class TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13> : ITypeMap
{
    public int Count => 14;


    public bool TryIndexOf<T>(out int index)
    {
        if (typeof(T) == typeof(T0))
        {
            index = 0;
            return true;
        }

        if (typeof(T) == typeof(T1))
        {
            index = 1;
            return true;
        }

        if (typeof(T) == typeof(T2))
        {
            index = 2;
            return true;
        }

        if (typeof(T) == typeof(T3))
        {
            index = 3;
            return true;
        }

        if (typeof(T) == typeof(T4))
        {
            index = 4;
            return true;
        }

        if (typeof(T) == typeof(T5))
        {
            index = 5;
            return true;
        }

        if (typeof(T) == typeof(T6))
        {
            index = 6;
            return true;
        }

        if (typeof(T) == typeof(T7))
        {
            index = 7;
            return true;
        }

        if (typeof(T) == typeof(T8))
        {
            index = 8;
            return true;
        }

        if (typeof(T) == typeof(T9))
        {
            index = 9;
            return true;
        }

        if (typeof(T) == typeof(T10))
        {
            index = 10;
            return true;
        }

        if (typeof(T) == typeof(T11))
        {
            index = 11;
            return true;
        }

        if (typeof(T) == typeof(T12))
        {
            index = 12;
            return true;
        }

        if (typeof(T) == typeof(T13))
        {
            index = 13;
            return true;
        }

        index = -1;
        return false;
    }
}

/// <summary>Type map for 15 types.</summary>
public sealed class TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14> : ITypeMap
{
    public int Count => 15;


    public bool TryIndexOf<T>(out int index)
    {
        if (typeof(T) == typeof(T0))
        {
            index = 0;
            return true;
        }

        if (typeof(T) == typeof(T1))
        {
            index = 1;
            return true;
        }

        if (typeof(T) == typeof(T2))
        {
            index = 2;
            return true;
        }

        if (typeof(T) == typeof(T3))
        {
            index = 3;
            return true;
        }

        if (typeof(T) == typeof(T4))
        {
            index = 4;
            return true;
        }

        if (typeof(T) == typeof(T5))
        {
            index = 5;
            return true;
        }

        if (typeof(T) == typeof(T6))
        {
            index = 6;
            return true;
        }

        if (typeof(T) == typeof(T7))
        {
            index = 7;
            return true;
        }

        if (typeof(T) == typeof(T8))
        {
            index = 8;
            return true;
        }

        if (typeof(T) == typeof(T9))
        {
            index = 9;
            return true;
        }

        if (typeof(T) == typeof(T10))
        {
            index = 10;
            return true;
        }

        if (typeof(T) == typeof(T11))
        {
            index = 11;
            return true;
        }

        if (typeof(T) == typeof(T12))
        {
            index = 12;
            return true;
        }

        if (typeof(T) == typeof(T13))
        {
            index = 13;
            return true;
        }

        if (typeof(T) == typeof(T14))
        {
            index = 14;
            return true;
        }

        index = -1;
        return false;
    }
}

/// <summary>Type map for 16 types.</summary>
public sealed class TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15> : ITypeMap
{
    public int Count => 16;


    public bool TryIndexOf<T>(out int index)
    {
        if (typeof(T) == typeof(T0))
        {
            index = 0;
            return true;
        }

        if (typeof(T) == typeof(T1))
        {
            index = 1;
            return true;
        }

        if (typeof(T) == typeof(T2))
        {
            index = 2;
            return true;
        }

        if (typeof(T) == typeof(T3))
        {
            index = 3;
            return true;
        }

        if (typeof(T) == typeof(T4))
        {
            index = 4;
            return true;
        }

        if (typeof(T) == typeof(T5))
        {
            index = 5;
            return true;
        }

        if (typeof(T) == typeof(T6))
        {
            index = 6;
            return true;
        }

        if (typeof(T) == typeof(T7))
        {
            index = 7;
            return true;
        }

        if (typeof(T) == typeof(T8))
        {
            index = 8;
            return true;
        }

        if (typeof(T) == typeof(T9))
        {
            index = 9;
            return true;
        }

        if (typeof(T) == typeof(T10))
        {
            index = 10;
            return true;
        }

        if (typeof(T) == typeof(T11))
        {
            index = 11;
            return true;
        }

        if (typeof(T) == typeof(T12))
        {
            index = 12;
            return true;
        }

        if (typeof(T) == typeof(T13))
        {
            index = 13;
            return true;
        }

        if (typeof(T) == typeof(T14))
        {
            index = 14;
            return true;
        }

        if (typeof(T) == typeof(T15))
        {
            index = 15;
            return true;
        }

        index = -1;
        return false;
    }
}

/// <summary>Type map for 17 types.</summary>
public sealed class TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16> : ITypeMap
{
    public int Count => 17;


    public bool TryIndexOf<T>(out int index)
    {
        if (typeof(T) == typeof(T0))
        {
            index = 0;
            return true;
        }

        if (typeof(T) == typeof(T1))
        {
            index = 1;
            return true;
        }

        if (typeof(T) == typeof(T2))
        {
            index = 2;
            return true;
        }

        if (typeof(T) == typeof(T3))
        {
            index = 3;
            return true;
        }

        if (typeof(T) == typeof(T4))
        {
            index = 4;
            return true;
        }

        if (typeof(T) == typeof(T5))
        {
            index = 5;
            return true;
        }

        if (typeof(T) == typeof(T6))
        {
            index = 6;
            return true;
        }

        if (typeof(T) == typeof(T7))
        {
            index = 7;
            return true;
        }

        if (typeof(T) == typeof(T8))
        {
            index = 8;
            return true;
        }

        if (typeof(T) == typeof(T9))
        {
            index = 9;
            return true;
        }

        if (typeof(T) == typeof(T10))
        {
            index = 10;
            return true;
        }

        if (typeof(T) == typeof(T11))
        {
            index = 11;
            return true;
        }

        if (typeof(T) == typeof(T12))
        {
            index = 12;
            return true;
        }

        if (typeof(T) == typeof(T13))
        {
            index = 13;
            return true;
        }

        if (typeof(T) == typeof(T14))
        {
            index = 14;
            return true;
        }

        if (typeof(T) == typeof(T15))
        {
            index = 15;
            return true;
        }

        if (typeof(T) == typeof(T16))
        {
            index = 16;
            return true;
        }

        index = -1;
        return false;
    }
}

/// <summary>Type map for 18 types.</summary>
public sealed class TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17> : ITypeMap
{
    public int Count => 18;


    public bool TryIndexOf<T>(out int index)
    {
        if (typeof(T) == typeof(T0))
        {
            index = 0;
            return true;
        }

        if (typeof(T) == typeof(T1))
        {
            index = 1;
            return true;
        }

        if (typeof(T) == typeof(T2))
        {
            index = 2;
            return true;
        }

        if (typeof(T) == typeof(T3))
        {
            index = 3;
            return true;
        }

        if (typeof(T) == typeof(T4))
        {
            index = 4;
            return true;
        }

        if (typeof(T) == typeof(T5))
        {
            index = 5;
            return true;
        }

        if (typeof(T) == typeof(T6))
        {
            index = 6;
            return true;
        }

        if (typeof(T) == typeof(T7))
        {
            index = 7;
            return true;
        }

        if (typeof(T) == typeof(T8))
        {
            index = 8;
            return true;
        }

        if (typeof(T) == typeof(T9))
        {
            index = 9;
            return true;
        }

        if (typeof(T) == typeof(T10))
        {
            index = 10;
            return true;
        }

        if (typeof(T) == typeof(T11))
        {
            index = 11;
            return true;
        }

        if (typeof(T) == typeof(T12))
        {
            index = 12;
            return true;
        }

        if (typeof(T) == typeof(T13))
        {
            index = 13;
            return true;
        }

        if (typeof(T) == typeof(T14))
        {
            index = 14;
            return true;
        }

        if (typeof(T) == typeof(T15))
        {
            index = 15;
            return true;
        }

        if (typeof(T) == typeof(T16))
        {
            index = 16;
            return true;
        }

        if (typeof(T) == typeof(T17))
        {
            index = 17;
            return true;
        }

        index = -1;
        return false;
    }
}

/// <summary>Type map for 19 types.</summary>
public sealed class
    TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18> : ITypeMap
{
    public int Count => 19;


    public bool TryIndexOf<T>(out int index)
    {
        if (typeof(T) == typeof(T0))
        {
            index = 0;
            return true;
        }

        if (typeof(T) == typeof(T1))
        {
            index = 1;
            return true;
        }

        if (typeof(T) == typeof(T2))
        {
            index = 2;
            return true;
        }

        if (typeof(T) == typeof(T3))
        {
            index = 3;
            return true;
        }

        if (typeof(T) == typeof(T4))
        {
            index = 4;
            return true;
        }

        if (typeof(T) == typeof(T5))
        {
            index = 5;
            return true;
        }

        if (typeof(T) == typeof(T6))
        {
            index = 6;
            return true;
        }

        if (typeof(T) == typeof(T7))
        {
            index = 7;
            return true;
        }

        if (typeof(T) == typeof(T8))
        {
            index = 8;
            return true;
        }

        if (typeof(T) == typeof(T9))
        {
            index = 9;
            return true;
        }

        if (typeof(T) == typeof(T10))
        {
            index = 10;
            return true;
        }

        if (typeof(T) == typeof(T11))
        {
            index = 11;
            return true;
        }

        if (typeof(T) == typeof(T12))
        {
            index = 12;
            return true;
        }

        if (typeof(T) == typeof(T13))
        {
            index = 13;
            return true;
        }

        if (typeof(T) == typeof(T14))
        {
            index = 14;
            return true;
        }

        if (typeof(T) == typeof(T15))
        {
            index = 15;
            return true;
        }

        if (typeof(T) == typeof(T16))
        {
            index = 16;
            return true;
        }

        if (typeof(T) == typeof(T17))
        {
            index = 17;
            return true;
        }

        if (typeof(T) == typeof(T18))
        {
            index = 18;
            return true;
        }

        index = -1;
        return false;
    }
}

/// <summary>Type map for 20 types.</summary>
public sealed class TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18,
    T19> : ITypeMap
{
    public int Count => 20;


    public bool TryIndexOf<T>(out int index)
    {
        if (typeof(T) == typeof(T0))
        {
            index = 0;
            return true;
        }

        if (typeof(T) == typeof(T1))
        {
            index = 1;
            return true;
        }

        if (typeof(T) == typeof(T2))
        {
            index = 2;
            return true;
        }

        if (typeof(T) == typeof(T3))
        {
            index = 3;
            return true;
        }

        if (typeof(T) == typeof(T4))
        {
            index = 4;
            return true;
        }

        if (typeof(T) == typeof(T5))
        {
            index = 5;
            return true;
        }

        if (typeof(T) == typeof(T6))
        {
            index = 6;
            return true;
        }

        if (typeof(T) == typeof(T7))
        {
            index = 7;
            return true;
        }

        if (typeof(T) == typeof(T8))
        {
            index = 8;
            return true;
        }

        if (typeof(T) == typeof(T9))
        {
            index = 9;
            return true;
        }

        if (typeof(T) == typeof(T10))
        {
            index = 10;
            return true;
        }

        if (typeof(T) == typeof(T11))
        {
            index = 11;
            return true;
        }

        if (typeof(T) == typeof(T12))
        {
            index = 12;
            return true;
        }

        if (typeof(T) == typeof(T13))
        {
            index = 13;
            return true;
        }

        if (typeof(T) == typeof(T14))
        {
            index = 14;
            return true;
        }

        if (typeof(T) == typeof(T15))
        {
            index = 15;
            return true;
        }

        if (typeof(T) == typeof(T16))
        {
            index = 16;
            return true;
        }

        if (typeof(T) == typeof(T17))
        {
            index = 17;
            return true;
        }

        if (typeof(T) == typeof(T18))
        {
            index = 18;
            return true;
        }

        if (typeof(T) == typeof(T19))
        {
            index = 19;
            return true;
        }

        index = -1;
        return false;
    }
}

/// <summary>Type map for 21 types.</summary>
public sealed class TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19,
    T20> : ITypeMap
{
    public int Count => 21;


    public bool TryIndexOf<T>(out int index)
    {
        if (typeof(T) == typeof(T0))
        {
            index = 0;
            return true;
        }

        if (typeof(T) == typeof(T1))
        {
            index = 1;
            return true;
        }

        if (typeof(T) == typeof(T2))
        {
            index = 2;
            return true;
        }

        if (typeof(T) == typeof(T3))
        {
            index = 3;
            return true;
        }

        if (typeof(T) == typeof(T4))
        {
            index = 4;
            return true;
        }

        if (typeof(T) == typeof(T5))
        {
            index = 5;
            return true;
        }

        if (typeof(T) == typeof(T6))
        {
            index = 6;
            return true;
        }

        if (typeof(T) == typeof(T7))
        {
            index = 7;
            return true;
        }

        if (typeof(T) == typeof(T8))
        {
            index = 8;
            return true;
        }

        if (typeof(T) == typeof(T9))
        {
            index = 9;
            return true;
        }

        if (typeof(T) == typeof(T10))
        {
            index = 10;
            return true;
        }

        if (typeof(T) == typeof(T11))
        {
            index = 11;
            return true;
        }

        if (typeof(T) == typeof(T12))
        {
            index = 12;
            return true;
        }

        if (typeof(T) == typeof(T13))
        {
            index = 13;
            return true;
        }

        if (typeof(T) == typeof(T14))
        {
            index = 14;
            return true;
        }

        if (typeof(T) == typeof(T15))
        {
            index = 15;
            return true;
        }

        if (typeof(T) == typeof(T16))
        {
            index = 16;
            return true;
        }

        if (typeof(T) == typeof(T17))
        {
            index = 17;
            return true;
        }

        if (typeof(T) == typeof(T18))
        {
            index = 18;
            return true;
        }

        if (typeof(T) == typeof(T19))
        {
            index = 19;
            return true;
        }

        if (typeof(T) == typeof(T20))
        {
            index = 20;
            return true;
        }

        index = -1;
        return false;
    }
}

/// <summary>Type map for 22 types.</summary>
public sealed class TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19,
    T20, T21> : ITypeMap
{
    public int Count => 22;


    public bool TryIndexOf<T>(out int index)
    {
        if (typeof(T) == typeof(T0))
        {
            index = 0;
            return true;
        }

        if (typeof(T) == typeof(T1))
        {
            index = 1;
            return true;
        }

        if (typeof(T) == typeof(T2))
        {
            index = 2;
            return true;
        }

        if (typeof(T) == typeof(T3))
        {
            index = 3;
            return true;
        }

        if (typeof(T) == typeof(T4))
        {
            index = 4;
            return true;
        }

        if (typeof(T) == typeof(T5))
        {
            index = 5;
            return true;
        }

        if (typeof(T) == typeof(T6))
        {
            index = 6;
            return true;
        }

        if (typeof(T) == typeof(T7))
        {
            index = 7;
            return true;
        }

        if (typeof(T) == typeof(T8))
        {
            index = 8;
            return true;
        }

        if (typeof(T) == typeof(T9))
        {
            index = 9;
            return true;
        }

        if (typeof(T) == typeof(T10))
        {
            index = 10;
            return true;
        }

        if (typeof(T) == typeof(T11))
        {
            index = 11;
            return true;
        }

        if (typeof(T) == typeof(T12))
        {
            index = 12;
            return true;
        }

        if (typeof(T) == typeof(T13))
        {
            index = 13;
            return true;
        }

        if (typeof(T) == typeof(T14))
        {
            index = 14;
            return true;
        }

        if (typeof(T) == typeof(T15))
        {
            index = 15;
            return true;
        }

        if (typeof(T) == typeof(T16))
        {
            index = 16;
            return true;
        }

        if (typeof(T) == typeof(T17))
        {
            index = 17;
            return true;
        }

        if (typeof(T) == typeof(T18))
        {
            index = 18;
            return true;
        }

        if (typeof(T) == typeof(T19))
        {
            index = 19;
            return true;
        }

        if (typeof(T) == typeof(T20))
        {
            index = 20;
            return true;
        }

        if (typeof(T) == typeof(T21))
        {
            index = 21;
            return true;
        }

        index = -1;
        return false;
    }
}

/// <summary>Type map for 23 types.</summary>
public sealed class TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19,
    T20, T21, T22> : ITypeMap
{
    public int Count => 23;


    public bool TryIndexOf<T>(out int index)
    {
        if (typeof(T) == typeof(T0))
        {
            index = 0;
            return true;
        }

        if (typeof(T) == typeof(T1))
        {
            index = 1;
            return true;
        }

        if (typeof(T) == typeof(T2))
        {
            index = 2;
            return true;
        }

        if (typeof(T) == typeof(T3))
        {
            index = 3;
            return true;
        }

        if (typeof(T) == typeof(T4))
        {
            index = 4;
            return true;
        }

        if (typeof(T) == typeof(T5))
        {
            index = 5;
            return true;
        }

        if (typeof(T) == typeof(T6))
        {
            index = 6;
            return true;
        }

        if (typeof(T) == typeof(T7))
        {
            index = 7;
            return true;
        }

        if (typeof(T) == typeof(T8))
        {
            index = 8;
            return true;
        }

        if (typeof(T) == typeof(T9))
        {
            index = 9;
            return true;
        }

        if (typeof(T) == typeof(T10))
        {
            index = 10;
            return true;
        }

        if (typeof(T) == typeof(T11))
        {
            index = 11;
            return true;
        }

        if (typeof(T) == typeof(T12))
        {
            index = 12;
            return true;
        }

        if (typeof(T) == typeof(T13))
        {
            index = 13;
            return true;
        }

        if (typeof(T) == typeof(T14))
        {
            index = 14;
            return true;
        }

        if (typeof(T) == typeof(T15))
        {
            index = 15;
            return true;
        }

        if (typeof(T) == typeof(T16))
        {
            index = 16;
            return true;
        }

        if (typeof(T) == typeof(T17))
        {
            index = 17;
            return true;
        }

        if (typeof(T) == typeof(T18))
        {
            index = 18;
            return true;
        }

        if (typeof(T) == typeof(T19))
        {
            index = 19;
            return true;
        }

        if (typeof(T) == typeof(T20))
        {
            index = 20;
            return true;
        }

        if (typeof(T) == typeof(T21))
        {
            index = 21;
            return true;
        }

        if (typeof(T) == typeof(T22))
        {
            index = 22;
            return true;
        }

        index = -1;
        return false;
    }
}

/// <summary>Type map for 24 types.</summary>
public sealed class TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19,
    T20, T21, T22, T23> : ITypeMap
{
    public int Count => 24;


    public bool TryIndexOf<T>(out int index)
    {
        if (typeof(T) == typeof(T0))
        {
            index = 0;
            return true;
        }

        if (typeof(T) == typeof(T1))
        {
            index = 1;
            return true;
        }

        if (typeof(T) == typeof(T2))
        {
            index = 2;
            return true;
        }

        if (typeof(T) == typeof(T3))
        {
            index = 3;
            return true;
        }

        if (typeof(T) == typeof(T4))
        {
            index = 4;
            return true;
        }

        if (typeof(T) == typeof(T5))
        {
            index = 5;
            return true;
        }

        if (typeof(T) == typeof(T6))
        {
            index = 6;
            return true;
        }

        if (typeof(T) == typeof(T7))
        {
            index = 7;
            return true;
        }

        if (typeof(T) == typeof(T8))
        {
            index = 8;
            return true;
        }

        if (typeof(T) == typeof(T9))
        {
            index = 9;
            return true;
        }

        if (typeof(T) == typeof(T10))
        {
            index = 10;
            return true;
        }

        if (typeof(T) == typeof(T11))
        {
            index = 11;
            return true;
        }

        if (typeof(T) == typeof(T12))
        {
            index = 12;
            return true;
        }

        if (typeof(T) == typeof(T13))
        {
            index = 13;
            return true;
        }

        if (typeof(T) == typeof(T14))
        {
            index = 14;
            return true;
        }

        if (typeof(T) == typeof(T15))
        {
            index = 15;
            return true;
        }

        if (typeof(T) == typeof(T16))
        {
            index = 16;
            return true;
        }

        if (typeof(T) == typeof(T17))
        {
            index = 17;
            return true;
        }

        if (typeof(T) == typeof(T18))
        {
            index = 18;
            return true;
        }

        if (typeof(T) == typeof(T19))
        {
            index = 19;
            return true;
        }

        if (typeof(T) == typeof(T20))
        {
            index = 20;
            return true;
        }

        if (typeof(T) == typeof(T21))
        {
            index = 21;
            return true;
        }

        if (typeof(T) == typeof(T22))
        {
            index = 22;
            return true;
        }

        if (typeof(T) == typeof(T23))
        {
            index = 23;
            return true;
        }

        index = -1;
        return false;
    }
}

/// <summary>Type map for 25 types.</summary>
public sealed class TypeMap<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19,
    T20, T21, T22, T23, T24> : ITypeMap
{
    public int Count => 25;


    public bool TryIndexOf<T>(out int index)
    {
        if (typeof(T) == typeof(T0))
        {
            index = 0;
            return true;
        }

        if (typeof(T) == typeof(T1))
        {
            index = 1;
            return true;
        }

        if (typeof(T) == typeof(T2))
        {
            index = 2;
            return true;
        }

        if (typeof(T) == typeof(T3))
        {
            index = 3;
            return true;
        }

        if (typeof(T) == typeof(T4))
        {
            index = 4;
            return true;
        }

        if (typeof(T) == typeof(T5))
        {
            index = 5;
            return true;
        }

        if (typeof(T) == typeof(T6))
        {
            index = 6;
            return true;
        }

        if (typeof(T) == typeof(T7))
        {
            index = 7;
            return true;
        }

        if (typeof(T) == typeof(T8))
        {
            index = 8;
            return true;
        }

        if (typeof(T) == typeof(T9))
        {
            index = 9;
            return true;
        }

        if (typeof(T) == typeof(T10))
        {
            index = 10;
            return true;
        }

        if (typeof(T) == typeof(T11))
        {
            index = 11;
            return true;
        }

        if (typeof(T) == typeof(T12))
        {
            index = 12;
            return true;
        }

        if (typeof(T) == typeof(T13))
        {
            index = 13;
            return true;
        }

        if (typeof(T) == typeof(T14))
        {
            index = 14;
            return true;
        }

        if (typeof(T) == typeof(T15))
        {
            index = 15;
            return true;
        }

        if (typeof(T) == typeof(T16))
        {
            index = 16;
            return true;
        }

        if (typeof(T) == typeof(T17))
        {
            index = 17;
            return true;
        }

        if (typeof(T) == typeof(T18))
        {
            index = 18;
            return true;
        }

        if (typeof(T) == typeof(T19))
        {
            index = 19;
            return true;
        }

        if (typeof(T) == typeof(T20))
        {
            index = 20;
            return true;
        }

        if (typeof(T) == typeof(T21))
        {
            index = 21;
            return true;
        }

        if (typeof(T) == typeof(T22))
        {
            index = 22;
            return true;
        }

        if (typeof(T) == typeof(T23))
        {
            index = 23;
            return true;
        }

        if (typeof(T) == typeof(T24))
        {
            index = 24;
            return true;
        }

        index = -1;
        return false;
    }
}

#endregion

#region IDataList - Type-erased list operations (replaces reflection)

/// <summary>
///     Interface for type-erased list operations, eliminating reflection needs.
/// </summary>
internal interface IDataList
{
    int Count { get; }
    int ElementSize { get; }

    // Issue #7 fix: Add reflection-free access for serialization performance
    Type ElementType { get; }
    void Clear();
    IDataList Clone();
    void ReorderByMapping(List<int> oldFromNew);
    void ReorderByInversePermutation(int[] inverse);
    object? GetItemAt(int index);
    void SetItemAt(int index, object? value);

    /// <summary>Adds an item to the list (type-erased).</summary>
    /// <param name="item">Item to add, must be compatible with the list's element type.</param>
    void AddItem(object? item);

    // Sub-entity extraction support: create empty list and copy items
    /// <summary>Creates an empty list of the same underlying type.</summary>
    IDataList CreateEmpty();

    /// <summary>Copies an item from source list at given index to this list.</summary>
    /// <param name="source">Source list (must be same underlying type).</param>
    /// <param name="index">Index in source list to copy from.</param>
    void AddFromSource(IDataList source, int index);

    /// <summary>Ensures the list has at least the specified capacity.</summary>
    void EnsureCapacity(int capacity);
}

/// <summary>
///     Generic implementation of IDataList for type-safe storage.
/// </summary>
internal sealed class DataList<T> : IDataList
{
    public DataList()
    {
        Items = new List<T>();
    }

    public DataList(int capacity)
    {
        Items = new List<T>(capacity);
    }

    public List<T> Items { get; }

    public T this[int index]
    {
        get => Items[index];
        set => Items[index] = value;
    }

    public int Count => Items.Count;

    public void AddItem(object? item)
    {
        if (item is T typedValue)
            Items.Add(typedValue);
        else if (item == null && (!typeof(T).IsValueType || Nullable.GetUnderlyingType(typeof(T)) != null))
            Items.Add(default!);
        else
            throw new ArgumentException(
                $"Cannot add {item?.GetType().Name ?? "null"} to DataList<{typeof(T).Name}>.");
    }

    public void Clear()
    {
        Items.Clear();
    }

    public IDataList Clone()
    {
        var clone = new DataList<T>(Items.Count);
        for (var i = 0; i < Items.Count; i++)
            clone.Items.Add(Items[i]);
        return clone;
    }

    public void ReorderByMapping(List<int> oldFromNew)
    {
        ArgumentNullException.ThrowIfNull(oldFromNew);

        if (oldFromNew.Count == 0 || Items.Count == 0) return;

        var newCount = Math.Min(oldFromNew.Count, Items.Count);

        // Validate all mapping indices are in valid range
        var span = CollectionsMarshal.AsSpan(oldFromNew);
        for (var newIdx = 0; newIdx < newCount; newIdx++)
        {
            var oldIdx = span[newIdx];
            if (oldIdx < 0 || oldIdx >= Items.Count)
                throw new ArgumentOutOfRangeException(
                    nameof(oldFromNew),
                    $"Mapping at index {newIdx} points to invalid source index {oldIdx}. " +
                    $"Valid range is [0, {Items.Count}).");
        }

        var temp = new T[newCount];

        for (var newIdx = 0; newIdx < newCount; newIdx++)
        {
            var oldIdx = span[newIdx];
            temp[newIdx] = Items[oldIdx];
        }

        Items.Clear();
        for (var i = 0; i < temp.Length; i++)
            Items.Add(temp[i]);
    }

    public void ReorderByInversePermutation(int[] inverse)
    {
        if (inverse.Length == 0 || Items.Count == 0) return;

        var count = Items.Count;
        if (inverse.Length < count)
            throw new ArgumentException(
                $"Inverse permutation length {inverse.Length} is less than item count {count}.",
                nameof(inverse));

        var temp = new T[count];

        for (var newIdx = 0; newIdx < count; newIdx++)
        {
            var oldIdx = inverse[newIdx];
            if (oldIdx < 0 || oldIdx >= count)
                throw new ArgumentOutOfRangeException(
                    nameof(inverse),
                    $"Inverse permutation at index {newIdx} points to invalid source index {oldIdx}. " +
                    $"Valid range is [0, {count}).");
            temp[newIdx] = Items[oldIdx];
        }

        for (var i = 0; i < count; i++)
            Items[i] = temp[i];
    }

    public int ElementSize
    {
        get
        {
            var type = typeof(T);
            if (type == typeof(int)) return 4;
            if (type == typeof(long)) return 8;
            if (type == typeof(float)) return 4;
            if (type == typeof(double)) return 8;
            if (type == typeof(bool)) return 1;
            if (type.IsValueType)
                try
                {
                    return Marshal.SizeOf(type);
                }
                catch
                {
                    return 16;
                }

            return IntPtr.Size;
        }
    }

    // Issue #7 fix: Reflection-free serialization support
    public Type ElementType => typeof(T);

    public object? GetItemAt(int index)
    {
        if (index < 0 || index >= Items.Count)
            throw new ArgumentOutOfRangeException(nameof(index));
        return Items[index];
    }

    public void SetItemAt(int index, object? value)
    {
        if (index < 0 || index >= Items.Count)
            throw new ArgumentOutOfRangeException(nameof(index));

        if (value == null)
        {
            if (!typeof(T).IsValueType || Nullable.GetUnderlyingType(typeof(T)) != null)
                Items[index] = default!;
            else
                throw new ArgumentNullException(nameof(value),
                    $"Cannot set null for non-nullable value type {typeof(T).Name}");
        }
        else if (value is T typedValue)
        {
            Items[index] = typedValue;
        }
        else
        {
            throw new ArgumentException(
                $"Value type {value.GetType().Name} is not compatible with {typeof(T).Name}",
                nameof(value));
        }
    }

    // Sub-entity extraction support
    public IDataList CreateEmpty()
    {
        return new DataList<T>();
    }

    public void AddFromSource(IDataList source, int index)
    {
        if (source is DataList<T> typedSource)
        {
            if (index < 0 || index >= typedSource.Items.Count)
                throw new ArgumentOutOfRangeException(nameof(index),
                    $"Index {index} is out of range for source list with {typedSource.Items.Count} items.");
            Items.Add(typedSource.Items[index]);
        }
        else
        {
            throw new InvalidOperationException(
                $"Source type mismatch: expected DataList<{typeof(T).Name}>, got {source.GetType().Name}");
        }
    }

    public void EnsureCapacity(int capacity)
    {
        if (Items.Capacity < capacity)
            Items.Capacity = capacity;
    }

    public void Add(T item)
    {
        Items.Add(item);
    }

    public Span<T> AsSpan()
    {
        return CollectionsMarshal.AsSpan(Items);
    }
}

#endregion

#region Utils

/// <summary>
///     Utility methods for list operations, set operations, and node mapping.
/// </summary>
public static class Utils
{
    #region Node Mapping

    /// <summary>
    ///     Validates that a kill list is properly normalized (sorted, no duplicates, in range).
    /// </summary>
    /// <param name="killList">The kill list to validate.</param>
    /// <param name="maxNodeValue">Maximum allowed node value.</param>
    /// <param name="throwOnError">If true, throws ArgumentException with detailed message on invalid input.</param>
    /// <returns>True if valid, false if invalid (when throwOnError is false).</returns>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool ValidateKillList(List<int> killList, int maxNodeValue, bool throwOnError = true)
    {
        if (killList.Count == 0)
            return true;

        var first = killList[0];
        if (first < 0 || first > maxNodeValue)
        {
            if (throwOnError)
                throw new ArgumentException(
                    $"Kill list contains out-of-range value {first} at index 0. " +
                    $"Valid range is [0, {maxNodeValue}].",
                    nameof(killList));
            return false;
        }

        for (var i = 1; i < killList.Count; i++)
        {
            var current = killList[i];
            var previous = killList[i - 1];

            if (current <= previous)
            {
                if (throwOnError)
                {
                    var msg = current == previous
                        ? $"Kill list contains duplicate value {current} at indices {i - 1} and {i}."
                        : $"Kill list is not sorted. Value {current} at index {i} is less than previous value {previous}.";
                    throw new ArgumentException(msg, nameof(killList));
                }

                return false;
            }

            if (current < 0 || current > maxNodeValue)
            {
                if (throwOnError)
                    throw new ArgumentException(
                        $"Kill list contains out-of-range value {current} at index {i}.",
                        nameof(killList));
                return false;
            }
        }

        return true;
    }

    public static (List<int> newNodesFromOld, List<int> oldNodesFromNew) GetNodeMapsFromKillList(
        int maxNodeValue, List<int> killList)
    {
        ArgumentNullException.ThrowIfNull(killList);
        if (maxNodeValue < 0) return (new List<int>(), new List<int>());

        // CRITICAL: Always validate kill list, not just in DEBUG (Review Issue #6)
        // Malformed kill lists cause silent corruption in release builds
        ValidateKillList(killList, maxNodeValue);

        var nodeCount = maxNodeValue + 1;
        var nodeExistsArray = ArrayPool<bool>.Shared.Rent(nodeCount);

        try
        {
            var nodeExists = nodeExistsArray.AsSpan(0, nodeCount);
            nodeExists.Fill(true);

            var killSpan = CollectionsMarshal.AsSpan(killList);
            for (var i = 0; i < killSpan.Length; i++)
            {
                var node = killSpan[i];
                if ((uint)node < (uint)nodeCount)
                    nodeExists[node] = false;
            }

            var oldToNew = new List<int>(nodeCount);
            var newToOld = new List<int>(nodeCount);
            var newIdx = 0;

            for (var oldIdx = 0; oldIdx < nodeCount; oldIdx++)
                if (nodeExists[oldIdx])
                {
                    oldToNew.Add(newIdx++);
                    newToOld.Add(oldIdx);
                }
                else
                {
                    oldToNew.Add(-1);
                }

            return (oldToNew, newToOld);
        }
        finally
        {
            ArrayPool<bool>.Shared.Return(nodeExistsArray);
        }
    }

    #endregion

    #region List Operations

    public static List<T> GetItemsAtIndices<T>(this List<T> list, List<int> indices)
    {
        ArgumentNullException.ThrowIfNull(list);
        ArgumentNullException.ThrowIfNull(indices);

        var result = new List<T>(indices.Count);
        var listSpan = CollectionsMarshal.AsSpan(list);
        var indicesSpan = CollectionsMarshal.AsSpan(indices);

        for (var i = 0; i < indicesSpan.Length; i++)
        {
            var index = indicesSpan[i];
            if ((uint)index >= (uint)listSpan.Length)
                throw new ArgumentOutOfRangeException(nameof(indices));
            result.Add(listSpan[index]);
        }

        return result;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void SortUnique<T>(this List<T> list) where T : IComparable<T>
    {
        ArgumentNullException.ThrowIfNull(list);
        if (list.Count <= 1) return;

        list.Sort();
        var span = CollectionsMarshal.AsSpan(list);
        var writeIndex = 1;

        for (var readIndex = 1; readIndex < span.Length; readIndex++)
            if (span[readIndex].CompareTo(span[writeIndex - 1]) != 0)
                span[writeIndex++] = span[readIndex];

        if (writeIndex < list.Count)
            list.RemoveRange(writeIndex, list.Count - writeIndex);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void InsertSorted(this List<int> sortedList, int value)
    {
        ArgumentNullException.ThrowIfNull(sortedList);

        var span = CollectionsMarshal.AsSpan(sortedList);
        int lo = 0, hi = span.Length;
        while (lo < hi)
        {
            var mid = lo + (hi - lo) / 2;
            if (span[mid] < value) lo = mid + 1;
            else hi = mid;
        }

        sortedList.Insert(lo, value);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static bool RemoveSorted(this List<int> sortedList, int value)
    {
        ArgumentNullException.ThrowIfNull(sortedList);

        var span = CollectionsMarshal.AsSpan(sortedList);
        int lo = 0, hi = span.Length - 1;
        while (lo <= hi)
        {
            var mid = lo + (hi - lo) / 2;
            if (span[mid] == value)
            {
                sortedList.RemoveAt(mid);
                return true;
            }

            if (span[mid] < value) lo = mid + 1;
            else hi = mid - 1;
        }

        return false;
    }

    #endregion

    #region Set Operations (Sorted Lists)

    [System.Diagnostics.Conditional("DEBUG")]
    private static void AssertSorted(List<int> list, string paramName)
    {
        var span = CollectionsMarshal.AsSpan(list);
        for (var i = 1; i < span.Length; i++)
            if (span[i] < span[i - 1])
                throw new ArgumentException(
                    $"List is not sorted: value {span[i]} at index {i} is less than {span[i - 1]} at index {i - 1}.",
                    paramName);
    }

    public static List<int> UnionSorted(List<int> a, List<int> b)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(b);
        AssertSorted(a, nameof(a));
        AssertSorted(b, nameof(b));

        var result = new List<int>(a.Count + b.Count);
        var aSpan = CollectionsMarshal.AsSpan(a);
        var bSpan = CollectionsMarshal.AsSpan(b);

        int i = 0, j = 0;
        while (i < aSpan.Length && j < bSpan.Length)
            if (aSpan[i] < bSpan[j])
            {
                result.Add(aSpan[i++]);
            }
            else if (aSpan[i] > bSpan[j])
            {
                result.Add(bSpan[j++]);
            }
            else
            {
                result.Add(aSpan[i++]);
                j++;
            }

        while (i < aSpan.Length) result.Add(aSpan[i++]);
        while (j < bSpan.Length) result.Add(bSpan[j++]);
        return result;
    }

    public static List<int> IntersectSorted(List<int> a, List<int> b)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(b);
        AssertSorted(a, nameof(a));
        AssertSorted(b, nameof(b));

        var result = new List<int>(Math.Min(a.Count, b.Count));
        var aSpan = CollectionsMarshal.AsSpan(a);
        var bSpan = CollectionsMarshal.AsSpan(b);

        int i = 0, j = 0;
        while (i < aSpan.Length && j < bSpan.Length)
            if (aSpan[i] < bSpan[j])
            {
                i++;
            }
            else if (aSpan[i] > bSpan[j])
            {
                j++;
            }
            else
            {
                result.Add(aSpan[i++]);
                j++;
            }

        return result;
    }

    public static List<int> DifferenceSorted(List<int> a, List<int> b)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(b);
        AssertSorted(a, nameof(a));
        AssertSorted(b, nameof(b));

        var result = new List<int>(a.Count);
        var aSpan = CollectionsMarshal.AsSpan(a);
        var bSpan = CollectionsMarshal.AsSpan(b);

        int i = 0, j = 0;
        while (i < aSpan.Length && j < bSpan.Length)
            if (aSpan[i] < bSpan[j])
            {
                result.Add(aSpan[i++]);
            }
            else if (aSpan[i] > bSpan[j])
            {
                j++;
            }
            else
            {
                i++;
                j++;
            }

        while (i < aSpan.Length) result.Add(aSpan[i++]);
        return result;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static bool ContainsSorted(List<int> sortedList, int value)
    {
        ArgumentNullException.ThrowIfNull(sortedList);
        if (sortedList.Count == 0) return false;

        var span = CollectionsMarshal.AsSpan(sortedList);
        int lo = 0, hi = span.Length - 1;
        while (lo <= hi)
        {
            var mid = lo + (hi - lo) / 2;
            if (span[mid] == value) return true;
            if (span[mid] < value) lo = mid + 1;
            else hi = mid - 1;
        }

        return false;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int BinarySearch(List<int> sortedList, int value)
    {
        ArgumentNullException.ThrowIfNull(sortedList);
        if (sortedList.Count == 0) return ~0;

        var span = CollectionsMarshal.AsSpan(sortedList);
        int lo = 0, hi = span.Length - 1;
        while (lo <= hi)
        {
            var mid = lo + (hi - lo) / 2;
            if (span[mid] == value) return mid;
            if (span[mid] < value) lo = mid + 1;
            else hi = mid - 1;
        }

        return ~lo;
    }

    #endregion

    #region Comparison

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int Compare<T>(List<T> first, List<T> second) where T : IComparable<T>
    {
        ArgumentNullException.ThrowIfNull(first);
        ArgumentNullException.ThrowIfNull(second);

        var firstSpan = CollectionsMarshal.AsSpan(first);
        var secondSpan = CollectionsMarshal.AsSpan(second);
        var minLength = Math.Min(firstSpan.Length, secondSpan.Length);

        for (var i = 0; i < minLength; i++)
        {
            var cmp = firstSpan[i].CompareTo(secondSpan[i]);
            if (cmp != 0) return cmp;
        }

        return first.Count.CompareTo(second.Count);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static bool AreEqual<T>(List<T> first, List<T> second) where T : IComparable<T>
    {
        ArgumentNullException.ThrowIfNull(first);
        ArgumentNullException.ThrowIfNull(second);
        if (first.Count != second.Count) return false;

        var firstSpan = CollectionsMarshal.AsSpan(first);
        var secondSpan = CollectionsMarshal.AsSpan(second);
        for (var i = 0; i < firstSpan.Length; i++)
            if (firstSpan[i].CompareTo(secondSpan[i]) != 0)
                return false;
        return true;
    }

    #endregion

    #region Min/Max/Sum

    public static int Min(List<int> list)
    {
        ArgumentNullException.ThrowIfNull(list);
        if (list.Count == 0) throw new InvalidOperationException("List is empty.");

        var span = CollectionsMarshal.AsSpan(list);
        var min = span[0];
        for (var i = 1; i < span.Length; i++)
            if (span[i] < min)
                min = span[i];
        return min;
    }

    public static int Max(List<int> list)
    {
        ArgumentNullException.ThrowIfNull(list);
        if (list.Count == 0) throw new InvalidOperationException("List is empty.");

        var span = CollectionsMarshal.AsSpan(list);
        var max = span[0];
        for (var i = 1; i < span.Length; i++)
            if (span[i] > max)
                max = span[i];
        return max;
    }

    public static long Sum(List<int> list)
    {
        ArgumentNullException.ThrowIfNull(list);
        var span = CollectionsMarshal.AsSpan(list);
        long sum = 0;
        for (var i = 0; i < span.Length; i++)
            sum += span[i];
        return sum;
    }

    public static int IndexOfMin(List<int> list)
    {
        ArgumentNullException.ThrowIfNull(list);
        if (list.Count == 0) throw new InvalidOperationException("List is empty.");

        var span = CollectionsMarshal.AsSpan(list);
        int minIdx = 0, min = span[0];
        for (var i = 1; i < span.Length; i++)
            if (span[i] < min)
            {
                min = span[i];
                minIdx = i;
            }

        return minIdx;
    }

    public static int IndexOfMax(List<int> list)
    {
        ArgumentNullException.ThrowIfNull(list);
        if (list.Count == 0) throw new InvalidOperationException("List is empty.");

        var span = CollectionsMarshal.AsSpan(list);
        int maxIdx = 0, max = span[0];
        for (var i = 1; i < span.Length; i++)
            if (span[i] > max)
            {
                max = span[i];
                maxIdx = i;
            }

        return maxIdx;
    }

    #endregion

    #region Copy Helpers (no LINQ)

    public static List<T> Copy<T>(List<T> source)
    {
        ArgumentNullException.ThrowIfNull(source);
        var result = new List<T>(source.Count);
        for (var i = 0; i < source.Count; i++)
            result.Add(source[i]);
        return result;
    }

    public static List<List<T>> DeepCopy<T>(List<List<T>> source)
    {
        ArgumentNullException.ThrowIfNull(source);
        var result = new List<List<T>>(source.Count);
        for (var i = 0; i < source.Count; i++)
        {
            var inner = new List<T>(source[i].Count);
            for (var j = 0; j < source[i].Count; j++)
                inner.Add(source[i][j]);
            result.Add(inner);
        }

        return result;
    }

    public static List<T> ToList<T>(HashSet<T> set)
    {
        ArgumentNullException.ThrowIfNull(set);
        var result = new List<T>(set.Count);
        foreach (var item in set)
            result.Add(item);
        return result;
    }

    public static List<T> ToSortedList<T>(HashSet<T> set) where T : IComparable<T>
    {
        var result = ToList(set);
        result.Sort();
        return result;
    }

    /// <summary>
    ///     Creates a list containing integers from 0 to count-1.
    /// </summary>
    public static List<int> Range(int count)
    {
        var result = new List<int>(count);
        for (var i = 0; i < count; i++)
            result.Add(i);
        return result;
    }

    #endregion

    #region Validation

    public static bool AreIndicesValid(List<int> indices, int collectionSize)
    {
        ArgumentNullException.ThrowIfNull(indices);
        var span = CollectionsMarshal.AsSpan(indices);
        for (var i = 0; i < span.Length; i++)
            if ((uint)span[i] >= (uint)collectionSize)
                return false;
        return true;
    }

    public static bool AreAllNonNegative(List<int> list)
    {
        ArgumentNullException.ThrowIfNull(list);
        var span = CollectionsMarshal.AsSpan(list);
        for (var i = 0; i < span.Length; i++)
            if (span[i] < 0)
                return false;
        return true;
    }

    #endregion
}

#endregion

#region Comparers

public sealed class ListComparer<T> : IComparer<List<T>> where T : IComparable<T>
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public int Compare(List<T>? x, List<T>? y)
    {
        if (ReferenceEquals(x, y)) return 0;
        if (x is null) return -1;
        if (y is null) return 1;
        return Utils.Compare(x, y);
    }
}

public sealed class ListEqualityComparer<T> : IEqualityComparer<List<T>> where T : IComparable<T>
{
    public bool Equals(List<T>? x, List<T>? y)
    {
        if (ReferenceEquals(x, y)) return true;
        if (x is null || y is null) return false;
        return Utils.AreEqual(x, y);
    }

    public int GetHashCode(List<T> obj)
    {
        ArgumentNullException.ThrowIfNull(obj);
        var hash = new HashCode();
        var span = CollectionsMarshal.AsSpan(obj);
        for (var i = 0; i < span.Length; i++)
            hash.Add(span[i]);
        return hash.ToHashCode();
    }
}

#endregion

#region ParallelConfig

/// <summary>
///     Global configuration for parallel operations.
///     Extended with GPU and MKL thread control.
/// </summary>
public static class ParallelConfig
{
    // ========================================================================
    // PRIVATE FIELDS
    // ========================================================================

    private static readonly object _lock = new();
    private static volatile int _maxDegreeOfParallelism = Environment.ProcessorCount;

    // GPU caching
    private static bool? _isGPUAvailable;
    private static volatile bool _skipGPUCheck; // Skip GPU check entirely if it causes crashes

    // MKL state
    private static int _mklNumThreads = Environment.ProcessorCount;
    private static IntPtr _mklHandle = IntPtr.Zero;
    private static MklSetNumThreadsDelegate? _mkl_set_num_threads;

    private static MklGetMaxThreadsDelegate? _mkl_get_max_threads;

    // P0-3 FIX: Split into two flags to allow retry after Cleanup()
    private static volatile bool _mklLoadAttempted;
    private static volatile bool _mklLoadedSuccessfully;

    // Debug logging

    // ========================================================================
    // STATIC CONSTRUCTOR
    // ========================================================================

    static ParallelConfig()
    {
        // Don't apply MKL threads on startup - wait until first use
        // This prevents segfaults if MKL is not available
        // ApplyMKLThreads() will be called on first MKL access
    }

    // ========================================================================
    // CPU PARALLELISM CONTROL
    // ========================================================================

    /// <summary>
    ///     Gets or sets whether debug output is enabled for diagnostics.
    ///     When true, diagnostic messages are written to System.Diagnostics.Debug.
    ///     Default: false
    /// </summary>
    public static bool EnableDebugOutput { get; set; }

    /// <summary>
    ///     Maximum degree of parallelism for all parallel operations.
    ///     Set to Environment.ProcessorCount by default.
    ///     Minimum value is 1 (no parallelism).
    /// </summary>
    public static int MaxDegreeOfParallelism
    {
        get => _maxDegreeOfParallelism;
        set
        {
            var newValue = Math.Max(1, value);
            if (newValue != _maxDegreeOfParallelism)
            {
                _maxDegreeOfParallelism = newValue;
                // P2.2 FIX: Invalidate cached ParallelOptions when config changes
                _cachedOptions = null;
            }
        }
    }

    // P2.2 FIX: Cache ParallelOptions to avoid allocation on every hot-path call
    private static volatile ParallelOptions? _cachedOptions;

    /// <summary>
    ///     Gets ParallelOptions configured with the global MaxDegreeOfParallelism.
    ///     P2.2 FIX: Cached to avoid allocation per call in hot paths.
    /// </summary>
    public static ParallelOptions Options
    {
        get
        {
            var opts = _cachedOptions;
            if (opts != null && opts.MaxDegreeOfParallelism == _maxDegreeOfParallelism)
                return opts;
            opts = new ParallelOptions { MaxDegreeOfParallelism = _maxDegreeOfParallelism };
            _cachedOptions = opts;
            return opts;
        }
    }

    // ========================================================================
    // GPU CONTROL
    // ========================================================================

    /// <summary>
    ///     Gets or sets whether to skip GPU availability checking entirely.
    ///     Set this to true if GPU detection causes crashes.
    ///     Default: false
    /// </summary>
    public static bool SkipGPUCheck
    {
        get => _skipGPUCheck;
        set
        {
            _skipGPUCheck = value;
            if (value)
                // Clear cached value so it won't be checked
                _isGPUAvailable = false;
            else
                // Allow re-checking after GPU checks are re-enabled
                _isGPUAvailable = null;
        }
    }

    /// <summary>
    ///     Gets or sets whether GPU acceleration is enabled globally.
    ///     When false, all operations use CPU regardless of GPU availability.
    ///     Default: true
    /// </summary>
    public static bool EnableGPU { get; set; } = true;

    /// <summary>
    ///     Returns true if GPU should be used (enabled AND available).
    /// </summary>
    public static bool UseGPU => EnableGPU && IsGPUAvailable;

    /// <summary>
    ///     Returns true if CUDA/cuSPARSE is available on this system.
    ///     Result is cached after first check.
    /// </summary>
    public static bool IsGPUAvailable
    {
        get
        {
            // If user explicitly disabled GPU checking, return false immediately
            if (_skipGPUCheck)
                return false;

            if (_isGPUAvailable.HasValue)
                return _isGPUAvailable.Value;

            lock (_lock)
            {
                // Re-check after acquiring lock
                if (_isGPUAvailable.HasValue)
                    return _isGPUAvailable.Value;

                var available = false;
                try
                {
                    // CRITICAL: Reflection can trigger static constructors that may crash
                    // Wrap the entire operation in multiple layers of protection
                    try
                    {
                        // P2-2 FIX: First try simple name (works if in same assembly)
                        var type = Type.GetType("Numerical.CuSparseBackend");

                        // P2-2 FIX: If not found, search all currently loaded assemblies
                        // This handles the case where CuSparseBackend is in a separate DLL
                        if (type == null)
                            try
                            {
                                foreach (var assembly in AppDomain.CurrentDomain.GetAssemblies())
                                {
                                    // Skip dynamic assemblies (can cause issues)
                                    if (assembly.IsDynamic)
                                        continue;

                                    type = assembly.GetType("Numerical.CuSparseBackend", false);
                                    if (type != null)
                                    {
                                        LogDebug($"Found CuSparseBackend in assembly: {assembly.FullName}");
                                        break;
                                    }
                                }
                            }
                            catch (Exception ex)
                            {
                                LogDebug($"Assembly search for GPU type failed: {ex.Message}");
                            }

                        if (type != null)
                            try
                            {
                                var prop = type.GetProperty("IsAvailable",
                                    BindingFlags.Public | BindingFlags.Static);
                                if (prop != null)
                                    try
                                    {
                                        var result = prop.GetValue(null);
                                        available = result is bool b && b;
                                    }
                                    catch (TargetInvocationException ex)
                                    {
                                        LogDebug(
                                            $"GPU property getter threw: {ex.InnerException?.Message ?? ex.Message}");
                                        available = false;
                                    }
                            }
                            catch (Exception ex)
                            {
                                LogDebug($"GPU property access failed: {ex.Message}");
                                available = false;
                            }
                    }
                    catch (TypeLoadException ex)
                    {
                        LogDebug($"GPU type load failed: {ex.Message}");
                        available = false;
                    }
                    catch (FileNotFoundException ex)
                    {
                        LogDebug($"GPU assembly not found: {ex.Message}");
                        available = false;
                    }
                }
                catch (Exception ex)
                {
                    LogDebug($"GPU availability check failed: {ex.Message}");
                    available = false;
                }

                _isGPUAvailable = available;
                return available;
            }
        }
    }

    // ========================================================================
    // MKL THREAD CONTROL
    // ========================================================================

    /// <summary>
    ///     Gets or sets the number of threads MKL uses internally.
    ///     Default: Environment.ProcessorCount
    ///     Valid range: 1 to 2*ProcessorCount (allows hyperthreading scenarios)
    /// </summary>
    public static int MKLNumThreads
    {
        get => _mklNumThreads;
        set
        {
            lock (_lock)
            {
                _mklNumThreads = Math.Clamp(value, 1, Environment.ProcessorCount * 2);
                ApplyMKLThreads();
            }
        }
    }

    /// <summary>
    ///     Gets the current number of threads MKL is using.
    ///     Returns null if MKL is not available.
    /// </summary>
    public static int? MKLCurrentThreads
    {
        get
        {
            EnsureMKLLoaded();

            // Only call if we have a valid function pointer
            if (_mkl_get_max_threads != null)
                try
                {
                    return _mkl_get_max_threads.Invoke();
                }
                catch (Exception ex)
                {
                    LogDebug($"Failed to get MKL threads: {ex.Message}");
                    return null;
                }

            return null;
        }
    }

    /// <summary>
    ///     Returns true if MKL is loaded and available.
    /// </summary>
    public static bool IsMKLAvailable => _mklLoadedSuccessfully && _mkl_set_num_threads != null;

    // ========================================================================
    // CONVENIENCE METHODS
    // ========================================================================

    /// <summary>
    ///     Shorthand for Environment.ProcessorCount.
    /// </summary>
    public static int ProcessorCount => Environment.ProcessorCount;

    /// <summary>
    ///     Writes a debug message if debug output is enabled.
    /// </summary>
    private static void LogDebug(string message)
    {
        if (EnableDebugOutput)
            Debug.WriteLine($"[ParallelConfig] {message}");
    }

    /// <summary>
    ///     Applies the current MKL thread count setting to the loaded MKL library.
    ///     Does nothing if MKL is not loaded or not available.
    /// </summary>
    private static void ApplyMKLThreads()
    {
        EnsureMKLLoaded();

        // Only call if we have a valid function pointer
        if (_mkl_set_num_threads != null)
            try
            {
                _mkl_set_num_threads.Invoke(_mklNumThreads);
            }
            catch (Exception ex)
            {
                LogDebug($"Failed to set MKL threads: {ex.Message}");
            }
    }

    /// <summary>
    ///     Ensures MKL library is loaded and function pointers are initialized.
    ///     Uses lazy loading with double-checked locking for thread safety.
    ///     Searches for platform-specific MKL library names and loads the first available.
    /// </summary>
    private static void EnsureMKLLoaded()
    {
        // P0-3 FIX: Use _mklLoadAttempted to allow retry after Cleanup()
        if (_mklLoadAttempted) return;

        lock (_lock)
        {
            if (_mklLoadAttempted) return;

            try
            {
                string[] mklNames;
                if (RuntimeInformation.IsOSPlatform(OSPlatform.Windows))
                    mklNames = new[] { "mkl_rt.2.dll", "mkl_rt.dll" };
                else if (RuntimeInformation.IsOSPlatform(OSPlatform.OSX))
                    mklNames = new[] { "libmkl_rt.2.dylib", "libmkl_rt.dylib" };
                else // Linux
                    mklNames = new[] { "libmkl_rt.so.2", "libmkl_rt.so" };

                var loadedSuccessfully = false;

                foreach (var name in mklNames)
                    try
                    {
                        if (NativeLibrary.TryLoad(name, out var handle))
                        {
                            // Free previous handle if we're retrying with a different library
                            if (_mklHandle != IntPtr.Zero && _mklHandle != handle)
                                try
                                {
                                    NativeLibrary.Free(_mklHandle);
                                    LogDebug("Freed previous MKL handle before loading new one");
                                }
                                catch (Exception ex)
                                {
                                    LogDebug($"Failed to free previous MKL handle: {ex.Message}");
                                }

                            _mklHandle = handle;

                            // Platform-specific function name priority (lowercase more common on Linux/macOS)
                            var setFuncNames = RuntimeInformation.IsOSPlatform(OSPlatform.Windows)
                                ? new[] { "MKL_Set_Num_Threads", "mkl_set_num_threads" }
                                : new[] { "mkl_set_num_threads", "MKL_Set_Num_Threads" };

                            var getFuncNames = RuntimeInformation.IsOSPlatform(OSPlatform.Windows)
                                ? new[] { "MKL_Get_Max_Threads", "mkl_get_max_threads" }
                                : new[] { "mkl_get_max_threads", "MKL_Get_Max_Threads" };

                            foreach (var funcName in setFuncNames)
                                if (NativeLibrary.TryGetExport(handle, funcName, out var setPtr))
                                    try
                                    {
                                        _mkl_set_num_threads =
                                            Marshal.GetDelegateForFunctionPointer<MklSetNumThreadsDelegate>(setPtr);
                                        break;
                                    }
                                    catch (Exception ex)
                                    {
                                        LogDebug($"Failed to create delegate for {funcName}: {ex.Message}");
                                    }

                            foreach (var funcName in getFuncNames)
                                if (NativeLibrary.TryGetExport(handle, funcName, out var getPtr))
                                    try
                                    {
                                        _mkl_get_max_threads =
                                            Marshal.GetDelegateForFunctionPointer<MklGetMaxThreadsDelegate>(getPtr);
                                        break;
                                    }
                                    catch (Exception ex)
                                    {
                                        LogDebug($"Failed to create delegate for {funcName}: {ex.Message}");
                                    }

                            if (_mkl_set_num_threads != null)
                            {
                                LogDebug($"MKL loaded successfully: {name}");
                                loadedSuccessfully = true;
                                break;
                            }

                            // Couldn't get function pointers, free this handle
                            try
                            {
                                NativeLibrary.Free(handle);
                                _mklHandle = IntPtr.Zero;
                            }
                            catch (Exception ex)
                            {
                                LogDebug($"Failed to free unusable MKL handle: {ex.Message}");
                            }
                        }
                    }
                    catch (Exception ex)
                    {
                        LogDebug($"Failed to load {name}: {ex.Message}");
                    }

                if (!loadedSuccessfully)
                {
                    LogDebug("MKL not found or could not be loaded");
                    _mklHandle = IntPtr.Zero;
                    _mkl_set_num_threads = null;
                    _mkl_get_max_threads = null;
                }

                // P0-3 FIX: Set both flags - attempted always true, successful only if loaded
                _mklLoadAttempted = true;
                _mklLoadedSuccessfully = loadedSuccessfully;
            }
            catch (Exception ex)
            {
                LogDebug($"MKL loading failed: {ex.GetType().Name}: {ex.Message}");
                _mklHandle = IntPtr.Zero;
                _mkl_set_num_threads = null;
                _mkl_get_max_threads = null;
                _mklLoadAttempted = true;
                _mklLoadedSuccessfully = false;
            }
        }
    }

    /// <summary>
    ///     Cleans up native MKL resources.
    ///     Call this before application exit if you need deterministic cleanup.
    ///     Note: This is typically not necessary as the OS will clean up on process exit.
    ///     After Cleanup(), MKL can be reloaded on next use (P0-3 fix).
    /// </summary>
    public static void Cleanup()
    {
        lock (_lock)
        {
            if (_mklHandle != IntPtr.Zero)
                try
                {
                    NativeLibrary.Free(_mklHandle);
                    LogDebug("MKL handle freed successfully");
                }
                catch (Exception ex)
                {
                    LogDebug($"MKL cleanup failed: {ex.GetType().Name}: {ex.Message}");
                }
                finally
                {
                    _mklHandle = IntPtr.Zero;
                    _mkl_set_num_threads = null;
                    _mkl_get_max_threads = null;
                    // P0-3 FIX: Reset load state to allow retry after cleanup
                    _mklLoadAttempted = false;
                    _mklLoadedSuccessfully = false;
                    LogDebug("MKL state reset - reload will be attempted on next use");
                }
        }
    }

    /// <summary>
    ///     Sets both CPU and MKL thread counts to the same value.
    /// </summary>
    /// <param name="numThreads">Number of threads (must be between 1 and 2*ProcessorCount)</param>
    /// <exception cref="ArgumentOutOfRangeException">Thrown if numThreads is out of valid range</exception>
    public static void SetAllThreads(int numThreads)
    {
        var maxAllowed = Environment.ProcessorCount * 2;
        if (numThreads < 1 || numThreads > maxAllowed)
            throw new ArgumentOutOfRangeException(nameof(numThreads),
                $"Thread count must be between 1 and {maxAllowed}, got {numThreads}");

        MaxDegreeOfParallelism = numThreads;
        MKLNumThreads = numThreads;
    }

    /// <summary>
    ///     Resets all settings to defaults.
    /// </summary>
    public static void Reset()
    {
        MaxDegreeOfParallelism = Environment.ProcessorCount;
        EnableGPU = true;
        EnableDebugOutput = false;
        _skipGPUCheck = false;
        _isGPUAvailable = null;
        MKLNumThreads = Environment.ProcessorCount;
    }

    /// <summary>
    ///     Gets a summary of current configuration.
    /// </summary>
    public static string GetSummary()
    {
        var gpuAvail = IsGPUAvailable; // Cache to avoid multiple reflection calls
        var mklThreads = MKLCurrentThreads;
        var mklStatus = IsMKLAvailable
            ? $"MKL={_mklNumThreads} (actual={mklThreads?.ToString() ?? "?"})"
            : "MKL=unavailable";

        return $"ParallelConfig: CPU={MaxDegreeOfParallelism}/{ProcessorCount}, " +
               $"GPU={EnableGPU} (Avail={gpuAvail}), " +
               $"{mklStatus}";
    }

    // ========================================================================
    // DELEGATES
    // ========================================================================

    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    private delegate void MklSetNumThreadsDelegate(int numThreads);

    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    private delegate int MklGetMaxThreadsDelegate();
}


#endregion
