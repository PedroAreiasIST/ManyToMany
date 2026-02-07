// EnsightWriter.cs - Export meshes to Ensight 6.0 format
// License: GPLv3

using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace Numerical;

/// <summary>
/// Exports SimplexMesh to Ensight 6.0 format (ASCII).
/// Supports both single-mesh and multi-mesh (time series) output.
/// </summary>
public static class EnsightWriter
{
    private static List<(string name, SimplexMesh mesh, double[,] coords, double[,]? displacement)> _meshCollection = new();
    
    /// <summary>
    /// Add mesh to collection for unified output.
    /// </summary>
    public static void AddMesh(string name, SimplexMesh mesh, double[,] coordinates, double[,]? displacement = null)
    {
        _meshCollection.Add((name, mesh, (double[,])coordinates.Clone(), 
            displacement != null ? (double[,])displacement.Clone() : null));
    }
    
    /// <summary>
    /// Write all collected meshes as a single Ensight case with multiple time steps.
    /// Each mesh appears as a separate time step for easy browsing in ParaView.
    /// </summary>
    public static void WriteAllMeshes(string basename)
    {
        if (_meshCollection.Count == 0)
        {
            Console.WriteLine("⚠️  No meshes collected for Ensight output");
            return;
        }
        
        string dir = Path.GetDirectoryName(basename);
        if (string.IsNullOrEmpty(dir))
            dir = ".";
        
        string name = Path.GetFileNameWithoutExtension(basename);
        if (string.IsNullOrEmpty(name))
            name = basename;
        
        Directory.CreateDirectory(dir);
        
        string caseFile = Path.Combine(dir, $"{name}.case");
        
        // Check if any mesh has displacement data
        bool hasDisplacement = _meshCollection.Any(m => m.displacement != null);
        
        // Write case file
        WriteCaseFileMultiStep(caseFile, name, _meshCollection.Count, hasDisplacement);
        
        // Write geometry files (one per time step)
        for (int i = 0; i < _meshCollection.Count; i++)
        {
            string geoFile = Path.Combine(dir, $"{name}_{i:D4}.geo");
            WriteGeometryFile(_meshCollection[i].mesh, _meshCollection[i].coords, 
                geoFile, _meshCollection[i].name);
            
            // Write displacement file if available
            if (_meshCollection[i].displacement != null)
            {
                string dispFile = Path.Combine(dir, $"{name}_{i:D4}.CrackOpening");
                WriteVectorFile(_meshCollection[i].mesh, _meshCollection[i].displacement!, 
                    dispFile, "CrackOpening");
            }
        }
        
        Console.WriteLine($"\n✓ Ensight output: {_meshCollection.Count} meshes → {name}.case");
        Console.WriteLine($"  Files created: {name}.case + {name}_XXXX.geo");
        if (hasDisplacement)
            Console.WriteLine($"  Displacement field: {name}_XXXX.CrackOpening (scale in Ensight!)");
        Console.WriteLine($"  Open in ParaView: File → Open → {name}.case");
        
        _meshCollection.Clear();
    }
    
    /// <summary>
    /// Save single mesh in Ensight 6.0 ASCII format.
    /// Creates: {basename}.case, {basename}.geo
    /// </summary>
    public static void SaveEnsight(SimplexMesh mesh, double[,] coordinates, string basename)
    {
        string dir = Path.GetDirectoryName(basename);
        if (string.IsNullOrEmpty(dir))
            dir = ".";
        
        string name = Path.GetFileNameWithoutExtension(basename);
        if (string.IsNullOrEmpty(name))
            name = basename;
        
        Directory.CreateDirectory(dir);
        
        string caseFile = Path.Combine(dir, $"{name}.case");
        string geoFile = Path.Combine(dir, $"{name}.geo");
        
        WriteCaseFileSingleStep(caseFile, name);
        WriteGeometryFile(mesh, coordinates, geoFile, name);
    }
    
    /// <summary>
    /// Write .case file for multi-step output.
    /// </summary>
    private static void WriteCaseFileMultiStep(string filename, string modelName, int numSteps, bool hasDisplacement = false)
    {
        using var writer = new StreamWriter(filename);
        
        writer.WriteLine("FORMAT");
        writer.WriteLine("type: ensight");  // Ensight 6.0
        writer.WriteLine();
        writer.WriteLine("GEOMETRY");
        writer.WriteLine($"model: 1 {numSteps} {modelName}_****.geo");
        writer.WriteLine();
        
        if (hasDisplacement)
        {
            writer.WriteLine("VARIABLE");
            writer.WriteLine($"vector per node: CrackOpening {modelName}_****.CrackOpening");
            writer.WriteLine();
        }
        
        writer.WriteLine("TIME");
        writer.WriteLine("time set: 1");
        writer.WriteLine($"number of steps: {numSteps}");
        writer.WriteLine("filename start number: 0");
        writer.WriteLine("filename increment: 1");
        writer.WriteLine("time values:");
        for (int i = 0; i < numSteps; i++)
        {
            writer.WriteLine($"{i,8}.0000");
        }
    }
    
    /// <summary>
    /// Write .case file for single-step output.
    /// </summary>
    private static void WriteCaseFileSingleStep(string filename, string modelName)
    {
        using var writer = new StreamWriter(filename);
        
        writer.WriteLine("FORMAT");
        writer.WriteLine("type: ensight");  // Ensight 6.0
        writer.WriteLine();
        writer.WriteLine("GEOMETRY");
        writer.WriteLine($"model: {modelName}.geo");
        writer.WriteLine();
        writer.WriteLine("TIME");
        writer.WriteLine("time set: 1");
        writer.WriteLine("number of steps: 1");
        writer.WriteLine("filename start number: 0");
        writer.WriteLine("filename increment: 1");
        writer.WriteLine("time values: 0.0");
    }
    
    /// <summary>
    /// Write .geo file (Ensight 6.0 geometry file).
    /// </summary>
    private static void WriteGeometryFile(SimplexMesh mesh, double[,] coordinates, string filename, string description)
    {
        using var writer = new StreamWriter(filename);
        var culture = CultureInfo.InvariantCulture;
        
        // Header
        writer.WriteLine($"{description}");
        writer.WriteLine("Generated by Numerical.SimplexMesh");
        writer.WriteLine("node id given");
        writer.WriteLine("element id given");
        
        // Part 1: All elements
        writer.WriteLine("part");
        writer.WriteLine("         1");
        writer.WriteLine($"{description}");
        writer.WriteLine("coordinates");
        
        // Node count
        int nNodes = mesh.Count<Node>();
        writer.WriteLine($"{nNodes,10}");
        
        // Node IDs
        for (int i = 0; i < nNodes; i++)
        {
            writer.WriteLine($"{i + 1,10}");
        }
        
        // X coordinates
        for (int i = 0; i < nNodes; i++)
        {
            writer.WriteLine(coordinates[i, 0].ToString("E5", culture).PadLeft(12));
        }
        
        // Y coordinates
        for (int i = 0; i < nNodes; i++)
        {
            writer.WriteLine(coordinates[i, 1].ToString("E5", culture).PadLeft(12));
        }
        
        // Z coordinates (all zero for 2D)
        for (int i = 0; i < nNodes; i++)
        {
            writer.WriteLine(0.0.ToString("E5", culture).PadLeft(12));
        }
        
        // Triangles (tria3)
        int nTri = mesh.Count<Tri3>();
        if (nTri > 0)
        {
            writer.WriteLine("tria3");
            writer.WriteLine($"{nTri,10}");
            
            // Element IDs
            for (int i = 0; i < nTri; i++)
            {
                writer.WriteLine($"{i + 1,10}");
            }
            
            // Connectivity (1-based indexing)
            for (int i = 0; i < nTri; i++)
            {
                var nodes = mesh.NodesOf<Tri3, Node>(i);
                writer.WriteLine($"{nodes[0] + 1,10}{nodes[1] + 1,10}{nodes[2] + 1,10}");
            }
        }
        
        // Quadrilaterals (quad4)
        int nQuad = mesh.Count<Quad4>();
        if (nQuad > 0)
        {
            writer.WriteLine("quad4");
            writer.WriteLine($"{nQuad,10}");
            
            // Element IDs (continue from triangles)
            for (int i = 0; i < nQuad; i++)
            {
                writer.WriteLine($"{nTri + i + 1,10}");
            }
            
            // Connectivity (1-based indexing)
            for (int i = 0; i < nQuad; i++)
            {
                var nodes = mesh.NodesOf<Quad4, Node>(i);
                writer.WriteLine($"{nodes[0] + 1,10}{nodes[1] + 1,10}{nodes[2] + 1,10}{nodes[3] + 1,10}");
            }
        }
    }
    
    /// <summary>
    /// Save mesh with scalar field data in Ensight format.
    /// Creates: {basename}.case, {basename}.geo, {basename}.{varname}
    /// </summary>
    public static void SaveEnsightWithScalar(SimplexMesh mesh, double[,] coordinates, 
        string basename, string variableName, double[] scalarData)
    {
        string dir = Path.GetDirectoryName(basename);
        if (string.IsNullOrEmpty(dir))
            dir = ".";
        
        string name = Path.GetFileNameWithoutExtension(basename);
        if (string.IsNullOrEmpty(name))
            name = basename;
        
        Directory.CreateDirectory(dir);
        
        string caseFile = Path.Combine(dir, $"{name}.case");
        string geoFile = Path.Combine(dir, $"{name}.geo");
        string varFile = Path.Combine(dir, $"{name}.{variableName}");
        
        WriteCaseFileWithVariable(caseFile, name, variableName);
        WriteGeometryFile(mesh, coordinates, geoFile, name);
        WriteScalarFile(mesh, scalarData, varFile, variableName);
    }
    
    private static void WriteCaseFileWithVariable(string filename, string modelName, string variableName)
    {
        using var writer = new StreamWriter(filename);
        
        writer.WriteLine("FORMAT");
        writer.WriteLine("type: ensight");  // Ensight 6.0
        writer.WriteLine();
        writer.WriteLine("GEOMETRY");
        writer.WriteLine($"model: {modelName}.geo");
        writer.WriteLine();
        writer.WriteLine("VARIABLE");
        writer.WriteLine($"scalar per node: {variableName} {modelName}.{variableName}");
        writer.WriteLine();
        writer.WriteLine("TIME");
        writer.WriteLine("time set: 1");
        writer.WriteLine("number of steps: 1");
        writer.WriteLine("filename start number: 0");
        writer.WriteLine("filename increment: 1");
        writer.WriteLine("time values: 0.0");
    }
    
    private static void WriteScalarFile(SimplexMesh mesh, double[] data, string filename, string description)
    {
        using var writer = new StreamWriter(filename);
        var culture = CultureInfo.InvariantCulture;
        
        writer.WriteLine($"{description}");
        writer.WriteLine("part");
        writer.WriteLine("         1");
        writer.WriteLine("coordinates");
        
        int nNodes = mesh.Count<Node>();
        for (int i = 0; i < nNodes; i++)
        {
            writer.WriteLine(data[i].ToString("E5", culture).PadLeft(12));
        }
    }
    
    /// <summary>
    /// Write vector field file (for displacement, velocity, etc.)
    /// </summary>
    private static void WriteVectorFile(SimplexMesh mesh, double[,] vectorData, string filename, string description)
    {
        using var writer = new StreamWriter(filename);
        
        writer.WriteLine($"{description}");
        writer.WriteLine("part");
        writer.WriteLine("         1");
        writer.WriteLine("coordinates");
        
        int nNodes = mesh.Count<Node>();
        
        // Write X components (exactly 12 chars: " 1.2345E+01" or "-1.2345E+01")
        for (int i = 0; i < nNodes; i++)
        {
            double val = vectorData[i, 0];
            string formatted = val.ToString("0.0000E+00", System.Globalization.CultureInfo.InvariantCulture);
            writer.WriteLine($"{formatted,12}");
        }
        
        // Write Y components
        for (int i = 0; i < nNodes; i++)
        {
            double val = vectorData[i, 1];
            string formatted = val.ToString("0.0000E+00", System.Globalization.CultureInfo.InvariantCulture);
            writer.WriteLine($"{formatted,12}");
        }
        
        // Write Z components
        for (int i = 0; i < nNodes; i++)
        {
            double val = vectorData[i, 2];
            string formatted = val.ToString("0.0000E+00", System.Globalization.CultureInfo.InvariantCulture);
            writer.WriteLine($"{formatted,12}");
        }
    }
}
