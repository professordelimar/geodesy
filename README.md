# Geodesy

Computational tools, numerical experiments and educational examples in **physical geodesy**, with emphasis on the Earth's gravity field, geopotential, geoid determination, height systems, reference surfaces and numerical methods.

## About

This repository brings together Python codes developed for computational applications in physical geodesy.

The main objective is to provide clear and reusable implementations for the mathematical and numerical treatment of the Earth's gravity field, reference ellipsoids, normal gravity, gravity anomalies and disturbances, geoid and quasigeoid determination, height systems, spherical harmonics and related geodetic computations.

The repository is intended for students, researchers and professionals interested in computational and physical geodesy.

## Scope

The repository may include implementations related to:

### Reference systems and reference surfaces

- Geodetic coordinate systems
- Cartesian and ellipsoidal coordinates
- Reference ellipsoids
- Geodetic latitude, longitude and ellipsoidal height
- Coordinate transformations
- Geocentric coordinates
- Local geodetic systems
- Reference frames
- Geodetic datums

### Gravity field of the Earth

- Newtonian potential
- Gravitational potential
- Centrifugal potential
- Gravity potential
- Gravity vector
- Level surfaces
- Plumb line
- Equipotential surfaces
- Geopotential numbers

### Normal gravity field

- Normal potential
- Reference ellipsoid
- Normal gravity
- Somigliana formula
- Normal gravity variation with latitude and height
- Clairaut-type relations

### Gravity anomalies and disturbances

- Gravity anomaly
- Gravity disturbance
- Free-air correction
- Bouguer correction
- Terrain effects
- Height anomaly
- Deflection of the vertical
- Gravity gradients

### Geoid and quasigeoid

- Geoid concept
- Quasigeoid concept
- Geoid undulation
- Height anomaly
- Bruns formula
- Stokes integral
- Molodensky theory
- Remove-compute-restore approaches
- Regional geoid modeling
- Numerical integration

### Height systems

- Ellipsoidal heights
- Orthometric heights
- Normal heights
- Dynamic heights
- Geopotential numbers
- Height transformations

### Spherical harmonics

- Legendre polynomials
- Associated Legendre functions
- Spherical harmonic expansion
- Fully normalized harmonics
- Global geopotential models
- Degree and order
- Spectral representation of the gravity field

### Numerical and computational methods

- Numerical integration
- Least-squares adjustment
- Regression
- Interpolation
- Gridding
- Fourier methods
- Fast Fourier Transform
- Spectral analysis
- Filtering
- Coordinate transformations
- Error propagation
- Numerical visualization

## Programming language

The primary programming language is **Python**.

Typical scientific and geodetic tools used in the repository include:

- NumPy
- SciPy
- Matplotlib
- Pandas
- SymPy
- Jupyter
- pyproj
- GeographicLib

Additional dependencies may be introduced when required by individual applications.

## Installation

Clone the repository:

```bash
git clone https://github.com/professordelimar/geodesy.git
```

Enter the repository:

```bash
cd geodesy
```

Install the basic dependencies:

```bash
pip install -r requirements.txt
```

## Running the codes

Python scripts can generally be executed with:

```bash
python filename.py
```

Jupyter notebooks can be opened with:

```bash
jupyter notebook
```

Specific instructions should be documented inside individual scripts or notebooks whenever additional data, parameters or dependencies are required.

## Scientific use

The codes in this repository are intended for computational experimentation, research, teaching and methodological development in physical geodesy.

Whenever a script reproduces, adapts or implements a published equation, geodetic model, transformation, gravity-field method or numerical procedure, the corresponding scientific reference should be identified in the source code, notebook or accompanying documentation.

## References

The conceptual basis of this repository is centered on standard works in physical geodesy, especially Torge and Müller, Heiskanen and Moritz, Hofmann-Wellenhof and Moritz, and related references.

See [REFERENCES.md](REFERENCES.md) for the main reference list.

## Citation

If this repository or its codes contribute to academic work, please cite the repository and the relevant scientific publications associated with the specific method being used.

Citation metadata are provided in [CITATION.cff](CITATION.cff).

## Related repositories

Computational methods in gravimetry and magnetometry are maintained separately in:

- https://github.com/professordelimar/gravmag

General applied geophysics methods are maintained in:

- https://github.com/professordelimar/geophysics

## Author

**Nelson Ribeiro-Filho**

## Copyright and ownership

Original source code, documentation and other original materials in this repository belong to:

**64.200.407 NELSON DE LIMA RIBEIRO FILHO - ME**  
**CNPJ 64.200.407/0001-45**

All rights reserved unless explicitly stated otherwise.

See [LICENSE](LICENSE) for details.

## Social links

- Instagram: https://www.instagram.com/professordelimar/
- YouTube: https://www.youtube.com/@professordelimar
