// Mm3Remeasure: the harness behind the revised ADES Tables 10-14 (see README.md).
//
//   dotnet <dll> [table10|table11|tables1213|table14|all]     (default: table14)

internal static class Program
{
    public static int Main(string[] args)
    {
        var which = args.Length > 0 ? args[0].ToLowerInvariant() : "table14";
        switch (which)
        {
            case "table10": Table10.Run(); break;
            case "table11": Table11.Run(); break;
            case "tables1213": Tables1213.Run(); break;
            case "table14": Table14.Run(); break;
            case "all":
                Table10.Run();
                Table11.Run();
                Tables1213.Run();
                Table14.Run();
                break;
            default:
                Console.Error.WriteLine($"unknown table '{which}'; expected table10|table11|tables1213|table14|all");
                return 1;
        }

        return 0;
    }
}
