package pl.poznan.put;

import org.apache.commons.lang3.tuple.Pair;
import org.immutables.value.Value;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

@Value.Immutable
public interface ClashGraph {
  private void searchDepthFirst(
      String vertex, Set<String> visited, Set<String> connectedComponent) {
    for (var next : adjacencyLists().getOrDefault(vertex, Collections.emptySet())) {
      if (!visited.contains(next)) {
        visited.add(next);
        connectedComponent.add(next);
        searchDepthFirst(next, visited, connectedComponent);
      }
    }
  }

  @Value.Lazy
  default Set<Pair<String, String>> adjacencyMatrix() {
    return adjacencyLists().entrySet().stream()
        .flatMap(entry -> entry.getValue().stream().map(s -> Pair.of(entry.getKey(), s)))
        .collect(Collectors.toSet());
  }

  Map<String, Set<String>> adjacencyLists();

  default Set<String> vertices() {
    return adjacencyLists().keySet();
  }

  default boolean isEmpty() {
    return adjacencyLists().isEmpty();
  }

  default boolean isClashFree() {
    return adjacencyLists().values().stream().map(Set::size).noneMatch(size -> size > 0);
  }

  default List<Set<String>> connectedComponents() {
    var visited = new HashSet<String>();
    var connectedComponents = new ArrayList<Set<String>>();

    for (var vertex : adjacencyLists().keySet()) {
      if (!visited.contains(vertex)) {
        visited.add(vertex);
        var connectedComponent = new HashSet<String>();
        connectedComponent.add(vertex);
        searchDepthFirst(vertex, visited, connectedComponent);
        connectedComponents.add(connectedComponent);
      }
    }

    return connectedComponents;
  }

  default ClashGraph induceSubgraph(Set<String> chains) {
    var map =
        adjacencyLists().entrySet().stream()
            .filter(entry -> chains.contains(entry.getKey()))
            .collect(
                Collectors.toMap(
                    Map.Entry::getKey,
                    entry ->
                        entry.getValue().stream()
                            .filter(chains::contains)
                            .collect(Collectors.toSet())));
    return ImmutableClashGraph.builder().adjacencyLists(map).build();
  }
}
