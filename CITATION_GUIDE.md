# Citation Guide

## How to Cite This Software

### Quick Citation (APA Style)

```
Ribeiro, D. (2025). High-Order Navier-Stokes Solvers: Newton-Raphson Finite 
Difference and Fourier Spectral Methods (Version 1.0.0) [Computer software]. 
https://github.com/diogoribeiro7/navier-stokes-solvers
```

### Quick Citation (BibTeX)

```bibtex
@software{ribeiro2025navierstokes,
  author = {Ribeiro, Diogo},
  title = {High-Order Navier-Stokes Solvers: Newton-Raphson Finite Difference and Fourier Spectral Methods},
  year = {2025},
  version = {1.0.0},
  url = {https://github.com/diogoribeiro7/navier-stokes-solvers},
  note = {Instituto Politécnico do Porto - ESMAD}
}
```

---

## About the Citation Files

### CITATION.cff

The `CITATION.cff` file follows the [Citation File Format](https://citation-file-format.github.io/) standard and enables:

- **GitHub Citation Widget**: Automatic "Cite this repository" button on GitHub
- **Zenodo Integration**: Automatic metadata import when archiving
- **Machine-Readable**: Tools can automatically extract citation information
- **Version Tracking**: Links citations to specific versions

**Features:**
- ✓ Author with ORCID
- ✓ Software metadata
- ✓ Key references
- ✓ DOI placeholder (update after Zenodo archiving)
- ✓ Multiple export formats supported

### .zenodo.json

The `.zenodo.json` file provides metadata for automatic Zenodo archiving:

- **Automatic Archiving**: GitHub-Zenodo integration uses this file
- **DOI Generation**: Each release gets a permanent DOI
- **Rich Metadata**: Comprehensive keywords, subjects, and references
- **Community Visibility**: Listed in relevant Zenodo communities

---

## Getting a DOI with Zenodo

### Step 1: Connect GitHub to Zenodo

1. Go to [Zenodo](https://zenodo.org/)
2. Log in with your GitHub account
3. Navigate to **Settings** → **GitHub**
4. Find `navier-stokes-solvers` and toggle it **ON**

### Step 2: Create a Release

```bash
# Tag the release
git tag -a v1.0.0 -m "Release v1.0.0: Initial public release"
git push origin v1.0.0

# Or create via GitHub web interface:
# Releases → Draft a new release → Create release
```

### Step 3: Update DOI

After Zenodo creates the DOI (usually within minutes):

1. Copy the DOI (format: `10.5281/zenodo.XXXXXXX`)
2. Update `CITATION.cff`:
   ```yaml
   identifiers:
     - type: doi
       value: 10.5281/zenodo.XXXXXXX  # Replace XXXXXXX
   ```
3. Update `.zenodo.json` if needed
4. Commit and push changes

### Step 4: Add DOI Badge to README

```markdown
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
```

---

## Citation Formats

### In Academic Papers

#### ACM Format
```
Diogo Ribeiro. 2025. High-Order Navier-Stokes Solvers: Newton-Raphson Finite 
Difference and Fourier Spectral Methods. Version 1.0.0. Retrieved from 
https://github.com/diogoribeiro7/navier-stokes-solvers
```

#### IEEE Format
```
D. Ribeiro, "High-Order Navier-Stokes Solvers: Newton-Raphson Finite Difference 
and Fourier Spectral Methods," 2025. [Online]. Available: 
https://github.com/diogoribeiro7/navier-stokes-solvers
```

#### Chicago Format
```
Ribeiro, Diogo. 2025. "High-Order Navier-Stokes Solvers: Newton-Raphson Finite 
Difference and Fourier Spectral Methods." Computer software. Version 1.0.0. 
https://github.com/diogoribeiro7/navier-stokes-solvers.
```

### BibTeX (Detailed)

```bibtex
@software{ribeiro2025navierstokes,
  author       = {Ribeiro, Diogo},
  title        = {{High-Order Navier-Stokes Solvers: Newton-Raphson 
                   Finite Difference and Fourier Spectral Methods}},
  month        = jan,
  year         = 2025,
  publisher    = {Zenodo},
  version      = {1.0.0},
  doi          = {10.5281/zenodo.XXXXXXX},
  url          = {https://doi.org/10.5281/zenodo.XXXXXXX},
  note         = {Instituto Politécnico do Porto - ESMAD},
  keywords     = {Navier-Stokes, CFD, spectral methods, finite difference}
}
```

### RIS Format

```ris
TY  - COMP
AU  - Ribeiro, Diogo
TI  - High-Order Navier-Stokes Solvers: Newton-Raphson Finite Difference and Fourier Spectral Methods
PY  - 2025
DA  - 2025/01/16
PB  - GitHub
UR  - https://github.com/diogoribeiro7/navier-stokes-solvers
DO  - 10.5281/zenodo.XXXXXXX
VL  - 1.0.0
KW  - Navier-Stokes equations
KW  - computational fluid dynamics
KW  - spectral methods
KW  - finite difference method
ER  -
```

### EndNote XML

```xml
<?xml version="1.0" encoding="UTF-8"?>
<xml>
  <records>
    <record>
      <ref-type name="Computer Program">9</ref-type>
      <contributors>
        <authors>
          <author>Ribeiro, Diogo</author>
        </authors>
      </contributors>
      <titles>
        <title>High-Order Navier-Stokes Solvers: Newton-Raphson Finite Difference and Fourier Spectral Methods</title>
      </titles>
      <dates>
        <year>2025</year>
      </dates>
      <urls>
        <related-urls>
          <url>https://github.com/diogoribeiro7/navier-stokes-solvers</url>
        </related-urls>
      </urls>
      <electronic-resource-num>10.5281/zenodo.XXXXXXX</electronic-resource-num>
    </record>
  </records>
</xml>
```

---

## Citing Specific Components

### If you use the Finite Difference Solver:

```
This work used the finite difference Navier-Stokes solver with Newton-Raphson 
iteration from Ribeiro (2025) [1], validated against the lid-driven cavity 
benchmarks of Ghia et al. (1982) [2].

[1] Ribeiro, D. (2025). High-Order Navier-Stokes Solvers...
[2] Ghia, U., et al. (1982). High-Re solutions for incompressible flow...
```

### If you use the Spectral Solver:

```
Simulations were performed using the Fourier spectral solver from Ribeiro (2025) [1], 
implementing the pseudospectral method with 2/3 dealiasing rule as described in 
Canuto et al. (2007) [2].

[1] Ribeiro, D. (2025). High-Order Navier-Stokes Solvers...
[2] Canuto, C., et al. (2007). Spectral Methods: Fundamentals...
```

### If you use the Test Suite:

```
Code validation followed the comprehensive test suite of Ribeiro (2025), including 
Taylor-Green vortex comparisons and conservation property verification.
```

---

## For Different Use Cases

### In a Research Paper (Methods Section)

```
The Navier-Stokes equations were solved using the open-source high-order solver 
of Ribeiro (2025), which implements both finite difference and spectral methods. 
We employed the spectral solver with 128×128 Fourier modes and 4th-order 
Runge-Kutta time integration, achieving spectral accuracy with L2 errors below 
10⁻⁶ when validated against analytical solutions.
```

### In a Thesis

```
Chapter 3: Numerical Methods

The numerical simulations in this thesis were performed using the Navier-Stokes 
solver developed by Ribeiro (2025) [ref], which provides two complementary 
approaches: a Newton-Raphson implicit finite difference scheme and a Fourier 
pseudospectral method. The software has been extensively validated against 
benchmark problems (Ghia et al., 1982) and analytical solutions (Taylor-Green 
vortex), with comprehensive test coverage exceeding 95%.
```

### In Documentation/Technical Reports

```
Software: High-Order Navier-Stokes Solvers v1.0.0 (Ribeiro, 2025)
Source: https://github.com/diogoribeiro7/navier-stokes-solvers
License: MIT
DOI: 10.5281/zenodo.XXXXXXX
```

### In Code/Software Dependencies

```python
# requirements.txt or similar
# Navier-Stokes solver - https://github.com/diogoribeiro7/navier-stokes-solvers
# Citation: Ribeiro, D. (2025). DOI: 10.5281/zenodo.XXXXXXX
```

---

## Verification of Citation Information

You can verify and export citation information using:

### GitHub Web Interface
1. Navigate to the repository
2. Click "Cite this repository" button (appears after adding CITATION.cff)
3. Choose your preferred format

### Command Line Tools

```bash
# Validate CITATION.cff
cffconvert --validate

# Convert to BibTeX
cffconvert -f bibtex -o citation.bib

# Convert to APA
cffconvert -f apalike

# Convert to EndNote
cffconvert -f endnote
```

Install `cffconvert`:
```bash
pip install cffconvert
```

---

## Updating Citation Information

### For New Versions

1. Update `version` in both `CITATION.cff` and `.zenodo.json`
2. Update `date-released` in `CITATION.cff`
3. Create new release on GitHub
4. Zenodo will create a version-specific DOI
5. Update README badges

### For New Contributors

Add to `CITATION.cff`:
```yaml
authors:
  - given-names: New
    family-names: Contributor
    email: contributor@example.com
    orcid: 'https://orcid.org/XXXX-XXXX-XXXX-XXXX'
```

---

## Questions?

- **GitHub Citation Format**: https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files
- **Zenodo Help**: https://help.zenodo.org/
- **ORCID**: https://orcid.org/

**Contact**: dfr@esmad.ipp.pt

---

**Note**: After receiving your DOI from Zenodo, remember to replace all instances of `XXXXXXX` in this guide and in the citation files!
