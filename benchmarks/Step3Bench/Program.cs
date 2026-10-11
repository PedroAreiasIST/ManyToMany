// Program.cs — modes of the step3 driver.
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
//   slab kind=... [options]      plane-strain slab examples (Slab.cs)
//   penny ...                    penny cracks in a cube: one inclined (Kassir-Sih) or two coplanar (coalescence)
//   surface ...                  semi-elliptical surface crack growing through the plate (breakthrough)
//   verify-count n seed cases    counting gates: random planar cracks, checker, Step2 oracle
//   criteria                     prints the kink angles of the three criteria for K_II/K_I ratios

using System.Globalization;

namespace Step3.Bench;

public static class Program
{
    public static int Main(string[] argv)
    {
        CultureInfo.DefaultThreadCurrentCulture = CultureInfo.InvariantCulture;
        if (argv.Length == 0)
        {
            Console.WriteLine("modes: slab, penny, surface, verify-count, criteria");
            return 1;
        }

        var a = new Args(argv.Skip(1));
        switch (argv[0])
        {
            case "slab":
                Slab.Run(a);
                return 0;
            case "criteria":
                CriteriaTable();
                return 0;
            case "verify-count":
                Verify.Count(a);
                return 0;
            case "penny":
                Penny.Run(a);
                return 0;
            case "surface":
                Surface.Run(a);
                return 0;
            default:
                Console.WriteLine($"unknown mode {argv[0]}");
                return 1;
        }
    }

    private static void CriteriaTable()
    {
        Console.WriteLine("beta  KII/KI   MTS      Schoellmann  Richard   (kink in degrees, K_III = 0, nu = 0.3)");
        foreach (var beta in new[] { 30, 40, 50, 60, 70, 80, 90 })
        {
            var (kI, kII) = Reference.InclinedCrack(1, 1, beta * Math.PI / 180);
            var mts = Criteria.MtsKink(kI, kII) * 180 / Math.PI;
            var (sch, _, _) = Criteria.Schoellmann(kI, kII, 0, 0.3);
            var (ric, _, _) = Criteria.Richard(kI, kII, 0);
            Console.WriteLine($"{beta,4}  {kII / kI,6:F3}  {mts,7:F2}  {sch * 180 / Math.PI,10:F2}  {ric * 180 / Math.PI,8:F2}");
        }
    }
}
