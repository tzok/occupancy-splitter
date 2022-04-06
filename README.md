# occupancy-splitter

## Requirements

- Java 11
- Maven

## Building

```sh
mvn package
```

## Usage

```sh
./occupancy-splitter -i <PATH-TO-MMCIF-FILE>
```

## Example

The example is based on structure with PDB id 488D described in:

> Capture and Visualization of a Catalytic RNA Enzyme-Product Complex Using Crystal Lattice Trapping and X-Ray Holographic Reconstruction. J.B. Murray et al. _Molecular Cell_. 2000. 5(2):279–287. doi:[10.1016/S1097-2765(00)80423-2](https://doi.org/10.1016/S1097-2765(00)80423-2)

![Input with overlaps](assets/488D.png)

```sh
$ ./occupancy-splitter -i 488D.cif
Output file: 488D-B-C.cif
Output file: 488D-D.cif
```

![Cleaved strands](assets/488D-B-C.png)

![Uncleaved substrate](assets/488D-D.png)
