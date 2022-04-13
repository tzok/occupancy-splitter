package pl.poznan.put;

import org.immutables.value.Value;

import java.util.Set;

@Value.Immutable
public interface ClashFreeChains {
  default int size() {
    return chains().size();
  }

  Set<String> chains();

  default boolean contains(String chain) {
    return chains().contains(chain);
  }

  default boolean containsAll(ClashFreeChains other) {
    return chains().containsAll(other.chains());
  }
}
