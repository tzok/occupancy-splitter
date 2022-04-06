package pl.poznan.put;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.immutables.value.Value;

@Value.Immutable
public interface Atom {
  String chain();

  int residue();

  String name();

  Vector3D coordinates();
}
