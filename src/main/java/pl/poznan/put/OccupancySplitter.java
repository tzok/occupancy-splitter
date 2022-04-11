package pl.poznan.put;

import com.google.common.collect.Sets;
import org.apache.commons.cli.*;
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
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class OccupancySplitter {
  private static final Logger LOGGER = LoggerFactory.getLogger(OccupancySplitter.class);
  private static final double PHOSPHORUS_VDW_RADIUS = 1.85;

  public static void main(final String[] args) throws IOException, ParseException {
    Locale.setDefault(Locale.US);

    // open the mmCIF file
    var commandLine = handleCommandLine(args);
    var inputPath = commandLine.getOptionValue("i");
    var mmCifFile = CifIO.readFromPath(Paths.get(inputPath)).as(StandardSchemata.MMCIF);
    var data = mmCifFile.getFirstBlock();

    // the main part
    var chainOccupancy = readOccupancyInfo(data);
    var chainsFractionalOccupancy = findChainsWithFractionalOccupancy(chainOccupancy);

    if (chainsFractionalOccupancy.size() < 2) {
      System.out.println("No clashes detected!");
      System.exit(0);
    }

    var chainAtoms = findPhosphorusAtoms(data, chainsFractionalOccupancy);
    var chainClashes = buildGraphOfClashes(chainAtoms);
    var allSolutions = findClashfreeChainCombinations(chainsFractionalOccupancy, chainClashes);
    var selectedSolutions = leaveOnlyLargestSubsets(allSolutions);
    createOutputFiles(
        inputPath, mmCifFile, chainOccupancy, chainsFractionalOccupancy, selectedSolutions);
  }

  private static void createOutputFiles(
      String inputPath,
      MmCifFile original,
      Map<String, Double> chainOccupancy,
      Set<String> chainsFractionalOccupancy,
      List<Set<String>> selectedSolutions)
      throws IOException {
    for (var acceptedChains : selectedSolutions) {
      var currentSolution =
          chainOccupancy.keySet().stream()
              .filter(
                  key -> acceptedChains.contains(key) || !chainsFractionalOccupancy.contains(key))
              .collect(Collectors.toSet());
      var outFile = OccupancySplitter.createCifCopyWithChainSelection(original, currentSolution);
      var base = inputPath.replace(".cif", "");
      var currentName = acceptedChains.stream().sorted().collect(Collectors.joining("-"));
      var outPath = Path.of(String.format("%s-%s.cif", base, currentName));
      CifIO.writeText(outFile, outPath);
      System.out.println("Output file: " + outPath);
    }
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

  private static List<Set<String>> leaveOnlyLargestSubsets(List<Set<String>> solutions) {
    // eliminate subsets which are already included in larger subsets
    // sort desc. by the subset size
    solutions.sort((t, t1) -> -Integer.compare(t.size(), t1.size()));

    List<Set<String>> accepted = new ArrayList<>();
    for (var candidate : solutions) {
      var isSuperset = accepted.stream().anyMatch(set -> set.containsAll(candidate));
      if (!isSuperset) {
        accepted.add(candidate);
      }
    }
    return accepted;
  }

  private static List<Set<String>> findClashfreeChainCombinations(
      Set<String> chainsFractionalOccupancy, Map<String, Set<String>> chainClashes) {
    // check if removing some chains with occupancy < 1.0 results in a clash-free structure
    // do it for all 2^n possible subsets of chains with occupancy < 1.0
    List<Set<String>> solutions = new ArrayList<>();
    for (var chainsToCheck : Sets.powerSet(chainsFractionalOccupancy)) {
      Map<String, Set<String>> copy =
          chainClashes.entrySet().stream()
              .collect(
                  Collectors.toMap(Map.Entry::getKey, entry -> new HashSet<>(entry.getValue())));

      for (var chain : chainsToCheck) {
        for (var clashed : copy.getOrDefault(chain, Collections.emptySet())) {
          copy.getOrDefault(clashed, Collections.emptySet()).remove(chain);
        }
        copy.remove(chain);
      }

      var clashesLeft = copy.values().stream().map(Set::size).anyMatch(i -> i > 0);

      if (!clashesLeft) {
        solutions.add(
            chainsFractionalOccupancy.stream()
                .filter(chain -> !chainsToCheck.contains(chain))
                .collect(Collectors.toSet()));
      }
    }
    return solutions;
  }

  private static Map<String, Set<String>> buildGraphOfClashes(
      Map<String, List<ImmutableAtom>> atoms) {
    return Sets.combinations(atoms.keySet(), 2).stream()
        .map(strings -> strings.toArray(new String[2]))
        .filter(chains -> OccupancySplitter.areClashing(atoms.get(chains[0]), atoms.get(chains[1])))
        .flatMap(strings -> Stream.of(strings, new String[] {strings[1], strings[0]}))
        .collect(
            Collectors.groupingBy(
                strings -> strings[0],
                Collectors.mapping(strings -> strings[1], Collectors.toSet())));
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

  private static Set<String> findChainsWithFractionalOccupancy(Map<String, Double> chainOccupancy) {
    return chainOccupancy.entrySet().stream()
        .filter(entry -> entry.getValue() < 1.0)
        .map(Map.Entry::getKey)
        .collect(Collectors.toSet());
  }

  private static Map<String, Double> readOccupancyInfo(MmCifBlock data) {
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

    var parser = new DefaultParser();
    var commandLine = parser.parse(options, args);

    if (!commandLine.hasOption("i")) {
      var helpFormatter = new HelpFormatter();
      helpFormatter.printHelp("occupancy-splitter", options);
      System.exit(1);
    }

    return commandLine;
  }
}
