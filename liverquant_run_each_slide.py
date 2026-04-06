# ============================================================
# LIVERQUANT - Run each .svs file one at a time
# Just change SLIDES_FOLDER to your folder path
# ============================================================

from liverquant import detect_fat_globules_wsi
from cv2geojson import export_annotations
from openslide import OpenSlide
import cv2 as cv
import numpy as np
import csv, os, glob

# ---- CHANGE THIS TO YOUR FOLDER ----
SLIDES_FOLDER = r"C:\Users\DEEPIKA\lipid quant"

# ---- OUTPUT (saved to Desktop automatically) ----
OUTPUT_CSV            = r"C:\Users\DEEPIKA\Desktop\PerinuclearLipid_Results.csv"
OUTPUT_GEOJSON_FOLDER = r"C:\Users\DEEPIKA\Desktop\geojson_outputs"

# ---- DETECTION SETTINGS ----
LOWERB              = [0, 0, 200]
UPPERB              = [180, 30, 255]
MIN_DIAMETER        = 8
MAX_DIAMETER        = 80
PERINUCLEAR_MAX_DIA = 25   # droplets <= 25 um = perinuclear (surrounds nucleus)
DOWNSAMPLE          = 2
CORES               = 4

# ============================================================

def process_one_slide(svs_path, writer, csv_file):
    """Processes a single .svs file and writes one row to CSV"""

    slide_name = os.path.basename(svs_path)
    print(f"\n--- Processing: {slide_name} ---")

    try:
        # Open slide
        slide = OpenSlide(svs_path)
        try:
            mpp = float(slide.properties.get('openslide.mpp-x', 0.5))
        except:
            mpp = 0.5
            print(f"  Warning: pixel size unknown, using {mpp} um/px")
        print(f"  Pixel size: {mpp} um/px")

        # Detect fat globules
        print("  Detecting fat globules...")
        spa, all_globules, roi, run_time = detect_fat_globules_wsi(
            slide,
            lowerb=LOWERB,
            upperb=UPPERB,
            tile_size=2048,
            overlap=128,
            downsample=DOWNSAMPLE,
            min_diameter=MIN_DIAMETER,
            max_diameter=MAX_DIAMETER,
            cores_num=CORES
        )
        total = len(all_globules)
        print(f"  Total globules detected: {total}")

        # Classify into perinuclear vs macrovesicular
        perinuclear, macrovesicular = [], []
        peri_area_px = 0
        tissue_area_px = sum(cv.contourArea(gc.contour) for gc in roi)

        for gc in all_globules:
            area_px   = cv.contourArea(gc.contour)
            radius_px = np.sqrt(area_px / np.pi)
            dia_um    = 2 * radius_px * mpp * DOWNSAMPLE
            if dia_um <= PERINUCLEAR_MAX_DIA:
                perinuclear.append(gc)
                peri_area_px += area_px
            else:
                macrovesicular.append(gc)

        # Calculate areas
        tissue_um2 = tissue_area_px * (mpp * DOWNSAMPLE) ** 2
        peri_um2   = peri_area_px   * (mpp * DOWNSAMPLE) ** 2
        peri_spa   = round((peri_um2 / tissue_um2 * 100) if tissue_um2 > 0 else 0, 2)
        peri_pct   = round((len(perinuclear) / total * 100) if total > 0 else 0, 1)

        # Print results
        print(f"  Perinuclear (<=25um): {len(perinuclear)}  ({peri_pct}%)")
        print(f"  Macrovesicular:       {len(macrovesicular)}")
        print(f"  Overall SPA:          {round(spa,2)}%")
        print(f"  Perinuclear SPA:      {peri_spa}%")
        print(f"  Run time:             {round(run_time,1)}s")

        # Save GeoJSON
        geojson_path = os.path.join(
            OUTPUT_GEOJSON_FOLDER,
            slide_name.replace('.svs', '_detections.geojson')
        )
        features = []
        for gc in roi:
            features.append(gc.export_feature(color=(0,0,255),   label='tissue'))
        for gc in perinuclear:
            features.append(gc.export_feature(color=(0,255,0),   label='perinuclear'))
        for gc in macrovesicular:
            features.append(gc.export_feature(color=(0,165,255), label='macrovesicular'))
        export_annotations(features, geojson_path)
        print(f"  GeoJSON saved: {os.path.basename(geojson_path)}")

        # Write to CSV
        writer.writerow([
            slide_name,
            total,
            len(perinuclear),
            len(macrovesicular),
            round(spa, 2),
            peri_spa,
            peri_pct,
            round(run_time, 1)
        ])
        csv_file.flush()
        print(f"  CSV updated.")

    except Exception as e:
        print(f"  ERROR: {e}")
        writer.writerow([slide_name, 'ERROR', '', '', '', '', '', ''])
        csv_file.flush()


# ============================================================
# MAIN - finds all .svs files and processes each one
# ============================================================

if __name__ == '__main__':

    # Find all .svs files
    svs_files = glob.glob(os.path.join(SLIDES_FOLDER, "*.svs"))

    if not svs_files:
        print("ERROR: No .svs files found in:", SLIDES_FOLDER)
        print("Check that SLIDES_FOLDER path is correct.")
    else:
        print(f"Found {len(svs_files)} slide(s):")
        for f in svs_files:
            print("  ", os.path.basename(f))

        # Create output folder
        os.makedirs(OUTPUT_GEOJSON_FOLDER, exist_ok=True)

        # Open CSV for writing
        write_header = not os.path.exists(OUTPUT_CSV)
        csv_file = open(OUTPUT_CSV, 'a', newline='')
        writer = csv.writer(csv_file)
        if write_header:
            writer.writerow([
                "Slide Name",
                "Total Fat Globules",
                "Perinuclear Globules",
                "Macrovesicular Globules",
                "Overall SPA (%)",
                "Perinuclear SPA (%)",
                "Perinuclear Globule %",
                "Run Time (s)"
            ])

        # --- LOOP: process each slide one by one ---
        for i, svs_path in enumerate(svs_files):
            print(f"\n[{i+1} of {len(svs_files)}]")
            process_one_slide(svs_path, writer, csv_file)

        csv_file.close()

        print("\n==========================================")
        print("ALL SLIDES DONE")
        print(f"Results CSV:  {OUTPUT_CSV}")
        print(f"GeoJSON files:{OUTPUT_GEOJSON_FOLDER}")
        print("==========================================")
