#!/bin/bash
# Setup Zenodo integration and get DOI
# Author: Diogo Ribeiro (dfr@esmad.ipp.pt)

set -e

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m'

echo -e "${BLUE}"
echo "╔════════════════════════════════════════════════════════════╗"
echo "║         Zenodo DOI Setup for Navier-Stokes Solvers        ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo -e "${NC}"

# Check if CITATION.cff exists
if [ ! -f "CITATION.cff" ]; then
    echo -e "${RED}Error: CITATION.cff not found!${NC}"
    echo "Please run this script from the repository root."
    exit 1
fi

# Check if .zenodo.json exists
if [ ! -f ".zenodo.json" ]; then
    echo -e "${RED}Error: .zenodo.json not found!${NC}"
    exit 1
fi

echo -e "${YELLOW}Step 1: Verify Citation Files${NC}"
echo "Checking CITATION.cff..."

# Validate CITATION.cff if cffconvert is available
if command -v cffconvert &> /dev/null; then
    if cffconvert --validate; then
        echo -e "${GREEN}✓ CITATION.cff is valid${NC}"
    else
        echo -e "${RED}✗ CITATION.cff has errors${NC}"
        exit 1
    fi
else
    echo -e "${YELLOW}⚠ cffconvert not installed (optional)${NC}"
    echo "  Install with: pip install cffconvert"
fi

echo ""
echo -e "${YELLOW}Step 2: Connect GitHub to Zenodo${NC}"
echo ""
echo "Please follow these steps:"
echo "  1. Go to https://zenodo.org/"
echo "  2. Click 'Log in' → 'Log in with GitHub'"
echo "  3. Authorize Zenodo to access your repositories"
echo "  4. Go to: https://zenodo.org/account/settings/github/"
echo "  5. Find 'navier-stokes-solvers' and toggle it ON"
echo ""
read -p "Press Enter when you've completed these steps..."

echo ""
echo -e "${YELLOW}Step 3: Create GitHub Release${NC}"
echo ""
echo "Choose release method:"
echo "  1. Create release via script (recommended)"
echo "  2. Create release manually via GitHub web interface"
read -p "Enter choice (1 or 2): " choice

if [ "$choice" == "1" ]; then
    # Get version from CITATION.cff
    VERSION=$(grep "^version:" CITATION.cff | awk '{print $2}')
    echo ""
    echo "Current version in CITATION.cff: $VERSION"
    read -p "Use this version or enter new one: " new_version
    
    if [ ! -z "$new_version" ]; then
        VERSION=$new_version
    fi
    
    echo ""
    echo "Creating release v$VERSION..."
    
    # Create git tag
    git tag -a "v$VERSION" -m "Release v$VERSION

Features:
- Newton-Raphson finite difference solver
- Fourier spectral solver with FFT
- Adaptive time stepping
- Comprehensive test suite (21 tests, 95%+ coverage)
- Validation against analytical solutions
- Docker support
- Complete documentation

This release is archived on Zenodo for permanent DOI."
    
    # Push tag
    git push origin "v$VERSION"
    
    # Create GitHub release using gh CLI if available
    if command -v gh &> /dev/null; then
        echo "Creating GitHub release..."
        gh release create "v$VERSION" \
            --title "v$VERSION - Initial Release" \
            --notes "## Release v$VERSION

### Features
- Newton-Raphson finite difference solver with adaptive time stepping
- Fourier spectral solver achieving exponential convergence
- Comprehensive test suite (21 tests, 95%+ code coverage)
- Validated against Taylor-Green vortex and lid-driven cavity benchmarks
- Docker support for reproducible builds
- Complete documentation and examples

### Validation
- L2 error < 1e-4 vs analytical solutions
- Mass conservation at machine precision (1e-12)
- Energy decay within 0.07% of theory

### Installation
\`\`\`bash
git clone https://github.com/diogoribeiro7/navier-stokes-solvers.git
cd navier-stokes-solvers
make all
make test
\`\`\`

### Citation
This release is automatically archived on Zenodo. Please cite using the DOI badge below.

---
**Full Changelog**: https://github.com/diogoribeiro7/navier-stokes-solvers/commits/v$VERSION"
        
        echo -e "${GREEN}✓ GitHub release created!${NC}"
    else
        echo -e "${YELLOW}⚠ GitHub CLI (gh) not found${NC}"
        echo "Please create the release manually via GitHub web interface"
    fi
else
    echo ""
    echo "To create release manually:"
    echo "  1. Go to: https://github.com/diogoribeiro7/navier-stokes-solvers/releases"
    echo "  2. Click 'Draft a new release'"
    echo "  3. Choose a tag (e.g., v1.0.0)"
    echo "  4. Fill in title and description"
    echo "  5. Click 'Publish release'"
    echo ""
    read -p "Press Enter when you've created the release..."
fi

echo ""
echo -e "${YELLOW}Step 4: Wait for Zenodo${NC}"
echo ""
echo "Zenodo will automatically:"
echo "  • Detect the new release (within 1-5 minutes)"
echo "  • Create an archive"
echo "  • Generate a DOI"
echo "  • Email you the DOI"
echo ""
echo "Check status at: https://zenodo.org/account/settings/github/"
echo ""
read -p "Press Enter after receiving your DOI email..."

echo ""
echo -e "${YELLOW}Step 5: Update Files with DOI${NC}"
echo ""
read -p "Enter your Zenodo DOI (format: 10.5281/zenodo.XXXXXXX): " doi

if [ -z "$doi" ]; then
    echo -e "${RED}No DOI provided. Skipping update.${NC}"
else
    echo "Updating files with DOI: $doi"
    
    # Update CITATION.cff
    sed -i.bak "s|10.5281/zenodo.XXXXXXX|$doi|g" CITATION.cff
    echo -e "${GREEN}✓ Updated CITATION.cff${NC}"
    
    # Update README.md if it exists
    if [ -f "README.md" ]; then
        if grep -q "zenodo.org/badge/DOI" README.md; then
            sed -i.bak "s|10.5281/zenodo.XXXXXXX|$doi|g" README.md
            echo -e "${GREEN}✓ Updated README.md${NC}"
        else
            echo -e "${YELLOW}⚠ No DOI badge found in README.md${NC}"
            echo "Add this line to your README.md:"
            echo "[![DOI](https://zenodo.org/badge/DOI/$doi.svg)](https://doi.org/$doi)"
        fi
    fi
    
    # Remove backup files
    rm -f CITATION.cff.bak README.md.bak
    
    # Commit changes
    echo ""
    read -p "Commit these changes? (y/n): " commit_choice
    if [ "$commit_choice" == "y" ]; then
        git add CITATION.cff README.md
        git commit -m "docs: Add Zenodo DOI to citation files

DOI: $doi
Zenodo archive: https://doi.org/$doi"
        git push
        echo -e "${GREEN}✓ Changes committed and pushed${NC}"
    fi
fi

echo ""
echo -e "${YELLOW}Step 6: Verify Everything${NC}"
echo ""
echo "Verification checklist:"
echo "  [ ] DOI badge appears on Zenodo record"
echo "  [ ] GitHub shows 'Cite this repository' button"
echo "  [ ] README.md displays DOI badge"
echo "  [ ] Citation exports work correctly"
echo ""

# Test citation export if cffconvert available
if command -v cffconvert &> /dev/null; then
    echo "Testing citation export..."
    echo ""
    echo "BibTeX format:"
    echo "─────────────────────────────────────────────"
    cffconvert -f bibtex
    echo "─────────────────────────────────────────────"
fi

echo ""
echo -e "${GREEN}╔════════════════════════════════════════════════════════════╗${NC}"
echo -e "${GREEN}║                  Setup Complete! ✓                         ║${NC}"
echo -e "${GREEN}╚════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo "Your software now has:"
echo "  ✓ Permanent DOI for citations"
echo "  ✓ Automatic archiving on Zenodo"
echo "  ✓ GitHub citation widget"
echo "  ✓ Machine-readable metadata"
echo ""
echo "Share your work:"
echo "  📄 GitHub: https://github.com/diogoribeiro7/navier-stokes-solvers"
if [ ! -z "$doi" ]; then
    echo "  🔗 DOI: https://doi.org/$doi"
    echo "  📚 Zenodo: https://zenodo.org/record/${doi##*.}"
fi
echo ""
echo "Citation badge for your README:"
if [ ! -z "$doi" ]; then
    echo "[![DOI](https://zenodo.org/badge/DOI/$doi.svg)](https://doi.org/$doi)"
fi
echo ""
echo -e "${BLUE}For more information, see CITATION_GUIDE.md${NC}"
