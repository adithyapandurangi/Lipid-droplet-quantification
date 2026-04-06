# Perinuclear Lipid Droplet Quantification
## H&E Rat Liver Whole Slide Image Analysis Pipeline

A two-tool pipeline for automated quantification of perinuclear lipid droplets in rat liver H&E stained whole slide images (.svs format).

---

## What This Does

- Detects fat globules (lipid vacuoles) in H&E liver sections
- Classifies droplets as **perinuclear** (small, ≤25 µm, surrounds nucleus) or **macrovesicular** (large, >25 µm, pushes nucleus to edge)
- Measures **Steatosis Proportionate Area (SPA%)** per slide
- Outputs per-nucleus perinuclear lipid measurements

---

## Pipeline Overview

```
.svs slides
    │
    ▼
Python / LiverQuant          ──►  PerinuclearLipid_Results.csv
(fat globule detection)      ──►  geojson_outputs/*.geojson
    │
    ▼
QuPath (groovy script)       ──►  QuPath_PerinuclearLipid_Results.csv
(nucleus detection +
 perinuclear matching)
```

---

## Files

| File | Description |
|---|---|
| `python/liverquant_run_each_slide.py` | Batch processes all .svs files, detects fat globules, exports GeoJSON + CSV |
| `qupath/QuPath_Fast_GeoJSON.groovy` | Loads GeoJSON from Python, detects nuclei, perinuclear spatial matching |
| `README_LipidQuantification.docx` | Full reproducible protocol with troubleshooting |

---

## Requirements

### Python
```bash
pip install "numpy==1.26.4" --force-reinstall
pip install liverquant openslide-bin openslide-python cv2geojson opencv-python
```

### QuPath
- QuPath v0.4.3+ (free at qupath.github.io)
- Set memory: Edit > Preferences > General > 24000 MB

---

## Quick Start

### Step 1 — Python
1. Open `liverquant_run_each_slide.py` in Spyder
2. Change `SLIDES_FOLDER` to your folder path
3. Press F5

### Step 2 — QuPath
1. Create a QuPath project and add all .svs files
2. Set each image to Brightfield (H&E)
3. Open `QuPath_Fast_GeoJSON.groovy` in script editor
4. Click Run > Run for Project

---

## Key Parameters

| Parameter | Default | Description |
|---|---|---|
| `PERINUCLEAR_MAX_DIA` | 25 µm | Droplets <= this = perinuclear |
| `DOWNSAMPLE` | 2 | 1=full res, 2=recommended, 4=testing |
| `EXPANSION_MICRONS` | 8.0 µm | Zone around nucleus (QuPath) |
| `MIN_DIAMETER` | 8 µm | Minimum droplet size |
| `MAX_DIAMETER` | 80 µm | Maximum droplet size |

---

## Output CSV Columns

**Python output:**
- `Total Fat Globules` — all detected lipid vacuoles
- `Perinuclear Globules` — small droplets surrounding nucleus
- `Overall SPA (%)` — % tissue area that is fat
- `Perinuclear SPA (%)` — % tissue area that is perinuclear fat

**QuPath output:**
- `% Nuclei With Lipid` — key steatosis measure
- `Avg Droplets per Nucleus` — lipid burden per hepatocyte
- `Avg Lipid Area per Nucleus (um2)` — area of lipid per cell

---

## Steatosis Grading

| SPA % | Grade | Meaning |
|---|---|---|
| < 5% | Grade 0 | Normal |
| 5–33% | Grade 1 | Mild steatosis |
| 33–66% | Grade 2 | Moderate steatosis |
| > 66% | Grade 3 | Severe steatosis |

---

## Troubleshooting

See `README_LipidQuantification.docx` for the full troubleshooting table covering all common errors.

---

## Citation / Acknowledgements

- Fat detection: [LiverQuant](https://github.com/mfarzi/liverquant) (Farzi et al.)
- Slide analysis: [QuPath](https://qupath.github.io) (Bankhead et al., 2017)
- Pipeline developed for rat liver steatosis quantification research

---

