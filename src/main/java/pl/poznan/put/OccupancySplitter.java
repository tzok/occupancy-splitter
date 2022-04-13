package pl.poznan.put;

import com.google.common.collect.Sets;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.lang3.NotImplementedException;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.rcsb.cif.CifBuilder;
import org.rcsb.cif.CifIO;
import org.rcsb.cif.model.FloatColumn;
import org.rcsb.cif.model.IntColumn;
import org.rcsb.cif.model.StrColumn;
import org.rcsb.cif.schema.StandardSchemata;
import org.rcsb.cif.schema.mm.MmCifBlock;
import org.rcsb.cif.schema.mm.MmCifFile;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class OccupancySplitter {
  private static final Logger LOGGER = LoggerFactory.getLogger(OccupancySplitter.class);
  private static final double PHOSPHORUS_VDW_RADIUS = 1.85;

  public static void main(String[] args) throws IOException, ParseException {
    Locale.setDefault(Locale.US);

    // open the mmCIF file
    var commandLine = handleCommandLine(args);
    var inputPath = commandLine.getOptionValue("i");
    var mmCifFile = CifIO.readFromPath(Paths.get(inputPath)).as(StandardSchemata.MMCIF);
    var data = mmCifFile.getFirstBlock();

    // the main part
    var chains = readChains(data);
    var chainsToCheck =
        commandLine.hasOption("all-chains")
            ? selectAllChains(chains)
            : selectChainsWithFractionalOccupancy(chains);
    var chainAtoms = findPhosphorusAtoms(data, chainsToCheck);
    var clashGraph = buildGraphOfClashes(chainAtoms);

    if (clashGraph.isClashFree()) {
      System.out.println("No clashes detected!");
      System.exit(0);
    }

    var solutionSets =
        clashGraph.connectedComponents().stream()
            .map(clashGraph::induceSubgraph)
            .map(OccupancySplitter::findClashFreeChainCombinations)
            .map(OccupancySplitter::leaveOnlyLargestSubsets)
            .collect(Collectors.toList());

    for (var clashFreeChains : unpackSolutionSets(solutionSets)) {
      createOutputFiles(inputPath, mmCifFile, chains, clashFreeChains);
    }
  }

  private static List<ClashFreeChains> unpackSolutionSets(
      List<List<ClashFreeChains>> solutionSets) {
    if (solutionSets.isEmpty()) {
      return Collections.emptyList();
    } else if (solutionSets.size() == 1) {
      return solutionSets.get(0);
    } else if (solutionSets.size() == 2) {
      List<ClashFreeChains> result = new ArrayList<>();
      for (var first : solutionSets.get(0)) {
        for (var second : solutionSets.get(1)) {
          result.add(
              ImmutableClashFreeChains.builder()
                  .addAllChains(first.chains())
                  .addAllChains(second.chains())
                  .build());
        }
      }
      return result;
    } else if (solutionSets.size() == 3) {
      List<ClashFreeChains> result = new ArrayList<>();
      for (var first : solutionSets.get(0)) {
        for (var second : solutionSets.get(1)) {
          for (var third : solutionSets.get(2)) {
            result.add(
                ImmutableClashFreeChains.builder()
                    .addAllChains(first.chains())
                    .addAllChains(second.chains())
                    .addAllChains(third.chains())
                    .build());
          }
        }
      }
      return result;
    }
    throw new NotImplementedException("Solution for 4+ connected components is not implemented");
  }

  private static Set<String> selectAllChains(Map<String, Double> chains) {
    return new HashSet<>(chains.keySet());
  }

  private static List<ClashFreeChains> leaveOnlyLargestSubsets(
      List<ClashFreeChains> clashFreeChains) {
    // eliminate subsets which are already included in larger subsets
    // sort desc. by the subset size
    clashFreeChains.sort((t, t1) -> -Integer.compare(t.size(), t1.size()));

    List<ClashFreeChains> accepted = new ArrayList<>();

    for (var candidate : clashFreeChains) {
      var isSuperset = accepted.stream().anyMatch(set -> set.containsAll(candidate));

      if (!isSuperset) {
        accepted.add(candidate);
      }
    }

    return accepted;
  }

  private static List<ClashFreeChains> findClashFreeChainCombinations(ClashGraph clashGraph) {
    return Sets.powerSet(clashGraph.vertices()).stream()
        .map(clashGraph::induceSubgraph)
        .filter(ClashGraph::isClashFree)
        .map(ClashGraph::vertices)
        .map(chains -> ImmutableClashFreeChains.builder().chains(chains).build())
        .collect(Collectors.toList());
  }

  private static ClashGraph buildGraphOfClashes(Map<String, List<ImmutableAtom>> atoms) {
    Map<String, Set<String>> map =
        atoms.size() >= 2
            ? Sets.combinations(atoms.keySet(), 2).stream()
                .map(strings -> strings.toArray(new String[2]))
                .filter(
                    chains ->
                        OccupancySplitter.areClashing(atoms.get(chains[0]), atoms.get(chains[1])))
                .flatMap(strings -> Stream.of(strings, new String[] {strings[1], strings[0]}))
                .collect(
                    Collectors.groupingBy(
                        strings -> strings[0],
                        Collectors.mapping(strings -> strings[1], Collectors.toSet())))
            : Collections.emptyMap();
    return ImmutableClashGraph.builder().adjacencyLists(map).build();
  }

  private static boolean areClashing(List<? extends Atom> chain1, List<? extends Atom> chain2) {
    for (Atom atom1 : chain1) {
      for (Atom atom2 : chain2) {
        if (atom1.coordinates().distanceSq(atom2.coordinates())
            < PHOSPHORUS_VDW_RADIUS * PHOSPHORUS_VDW_RADIUS) {
          return true;
        }
      }
    }
    return false;
  }

  private static Map<String, List<ImmutableAtom>> findPhosphorusAtoms(
      MmCifBlock data, Set<String> chainsFractionalOccupancy) {
    var atomSite = data.getAtomSite();
    var labelAsymId = atomSite.getLabelAsymId();
    var labelSeqId = atomSite.getLabelSeqId();
    var labelAtomId = atomSite.getLabelAtomId();
    var cartnX = atomSite.getCartnX();
    var cartnY = atomSite.getCartnY();
    var cartnZ = atomSite.getCartnZ();

    // phosphorus atom coordinates for every residue in every chain
    return IntStream.range(0, atomSite.getRowCount())
        .filter(i -> chainsFractionalOccupancy.contains(labelAsymId.get(i)))
        .filter(i -> "P".equals(labelAtomId.get(i)))
        .mapToObj(
            i ->
                ImmutableAtom.builder()
                    .chain(labelAsymId.get(i))
                    .residue(labelSeqId.get(i))
                    .name(labelAtomId.get(i))
                    .coordinates(new Vector3D(cartnX.get(i), cartnY.get(i), cartnZ.get(i)))
                    .build())
        .collect(Collectors.groupingBy(ImmutableAtom::chain));
  }

  private static Set<String> selectChainsWithFractionalOccupancy(Map<String, Double> chains) {
    return chains.entrySet().stream()
        .filter(entry -> entry.getValue() < 1.0)
        .map(Map.Entry::getKey)
        .collect(Collectors.toSet());
  }

  private static Map<String, Double> readChains(MmCifBlock data) {
    var atomSite = data.getAtomSite();
    var labelAsymId = atomSite.getLabelAsymId();
    var occupancy = atomSite.getOccupancy();

    return IntStream.range(0, atomSite.getRowCount())
        .mapToObj(i -> Pair.of(labelAsymId.get(i), occupancy.get(i)))
        .distinct()
        .collect(
            Collectors.toMap(
                Pair::getKey,
                Pair::getValue,
                (occupancy1, occupancy2) -> occupancy1 < occupancy2 ? occupancy1 : occupancy2));
  }

  private static CommandLine handleCommandLine(String[] args) throws ParseException {
    var options = new Options();
    options.addOption("i", "input", true, "(required) path to input file");
    options.addOption(
        "a",
        "all-chains",
        false,
        "(optional) check all chains for clashes, not only those with occupancy < 1.0");

    var parser = new DefaultParser();
    var commandLine = parser.parse(options, args);

    if (!commandLine.hasOption("i")) {
      var helpFormatter = new HelpFormatter();
      helpFormatter.printHelp("occupancy-splitter", options);
      System.exit(1);
    }

    return commandLine;
  }

  private static void createOutputFiles(
      String inputPath,
      MmCifFile original,
      Map<String, Double> chains,
      ClashFreeChains clashFreeChains)
      throws IOException {
    var currentSolution =
        chains.keySet().stream()
            .filter(key -> clashFreeChains.contains(key) || chains.get(key) == 1.0)
            .collect(Collectors.toSet());
    var outFile = OccupancySplitter.createCifCopyWithChainSelection(original, currentSolution);
    var base = inputPath.replace(".cif", "");
    var currentName = clashFreeChains.chains().stream().sorted().collect(Collectors.joining("-"));
    var outPath = Path.of(String.format("%s-%s.cif", base, currentName));
    CifIO.writeText(outFile, outPath);
    System.out.println("Output file: " + outPath);
  }

  private static MmCifFile createCifCopyWithChainSelection(
      MmCifFile original, Set<String> acceptedChains) {
    // gather indices of atoms which should stay
    var originalAtomSite = original.getFirstBlock().getAtomSite();
    var acceptedIndices =
        IntStream.range(0, originalAtomSite.getRowCount())
            .filter(i -> acceptedChains.contains(originalAtomSite.getLabelAsymId().get(i)))
            .boxed()
            .collect(Collectors.toList());

    var fileBuilder = CifBuilder.enterFile(StandardSchemata.MMCIF);

    for (var block : original.getBlocks()) {
      var blockBuilder = fileBuilder.enterBlock(block.getBlockHeader());

      for (var categoryEntry : block.getCategories().entrySet()) {
        var categoryName = categoryEntry.getKey();
        var category = categoryEntry.getValue();
        var categoryBuilder = blockBuilder.enterCategory(categoryName);

        if ("atom_site".equals(categoryName)) {
          // in atom_site, copy only selected lines
          for (var columnEntry : category.getColumns().entrySet()) {
            var columnName = columnEntry.getKey();
            var column = columnEntry.getValue();

            if (column instanceof StrColumn) {
              categoryBuilder
                  .enterStrColumn(columnName)
                  .add(
                      acceptedIndices.stream()
                          .map(i -> ((StrColumn) column).get(i))
                          .toArray(String[]::new))
                  .leaveColumn();
            } else if (column instanceof IntColumn) {
              categoryBuilder
                  .enterIntColumn(columnName)
                  .add(
                      acceptedIndices.stream().mapToInt(i -> ((IntColumn) column).get(i)).toArray())
                  .leaveColumn();
            } else if (column instanceof FloatColumn) {
              categoryBuilder
                  .enterFloatColumn(columnName)
                  .add(
                      acceptedIndices.stream()
                          .mapToDouble(i -> ((FloatColumn) column).get(i))
                          .toArray())
                  .leaveColumn();
            } else {
              OccupancySplitter.LOGGER.error("Invalid column: {}", column);
            }
          }
        } else {
          // in categories other than atom_site, copy everything as it is
          for (var columnEntry : category.getColumns().entrySet()) {
            var columnName = columnEntry.getKey();
            var column = columnEntry.getValue();

            if (column instanceof StrColumn) {
              categoryBuilder
                  .enterStrColumn(columnName)
                  .add(((StrColumn) column).values().toArray(String[]::new))
                  .leaveColumn();
            } else if (column instanceof IntColumn) {
              categoryBuilder
                  .enterIntColumn(columnName)
                  .add(((IntColumn) column).values().toArray())
                  .leaveColumn();
            } else if (column instanceof FloatColumn) {
              categoryBuilder
                  .enterFloatColumn(columnName)
                  .add(((FloatColumn) column).values().toArray())
                  .leaveColumn();
            } else {
              OccupancySplitter.LOGGER.error("Invalid column: {}", column);
            }
          }
        }

        categoryBuilder.leaveCategory();
      }

      blockBuilder.leaveBlock();
    }

    return fileBuilder.leaveFile();
  }
}
