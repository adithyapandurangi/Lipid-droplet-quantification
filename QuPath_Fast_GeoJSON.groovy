// ============================================================
// QUPATH - Fast Perinuclear Lipid using Python GeoJSON
// Loads fat globules from Python/LiverQuant GeoJSON files
// Uses spatial index for fast perinuclear matching
// ============================================================
// HOW TO USE:
//   1. Make sure Python script has finished and GeoJSON files
//      are in Desktop/geojson_outputs/
//   2. Open QuPath project with all .svs files
//   3. Set each image to Brightfield (H&E)
//   4. Automate > Script Editor > File > Open this script
//   5. Run > Run for Project
// ============================================================

def imageData = getCurrentImageData()
def server    = imageData.getServer()
def cal       = server.getPixelCalibration()
def pixelSize = cal.getAveragedPixelSizeMicrons()
def slideName = server.getMetadata().getName()
def outputCSV = System.getProperty("user.home") + "/Desktop/QuPath_PerinuclearLipid_Results.csv"

// GeoJSON folder — where Python saved its output
def GEOJSON_FOLDER = System.getProperty("user.home") + "/Desktop/geojson_outputs"

// Perinuclear expansion zone in microns
double EXPANSION_MICRONS = 8.0

println "\n=============================="
println "Processing: " + slideName
println "=============================="

// STEP 1: Detect tissue
println "Step 1: Detecting tissue..."
clearAllObjects()
runPlugin('qupath.imagej.detect.tissue.SimpleTissueDetection2',
    '{"threshold": 230, "requestedPixelSizeMicrons": 5.0,' +
    '"minAreaMicrons": 10000.0, "maxHoleAreaMicrons": 1000000.0,' +
    '"darkBackground": false, "smoothImage": true,' +
    '"medianCleanup": true, "dilateBoundaries": false,' +
    '"smoothCoordinates": true, "excludeOnBoundary": false,' +
    '"singleAnnotation": true}')

def tissueAnns = getAnnotationObjects()
if (tissueAnns.isEmpty()) {
    println "ERROR: No tissue detected. Skipping."
    return
}
double tissuePx = 0
tissueAnns.each { tissuePx += it.getROI().getArea() }
double tissueArea = (!Double.isNaN(pixelSize) && pixelSize > 0) ?
    tissuePx * pixelSize * pixelSize : tissuePx
println "  Tissue area: " + String.format("%.1f", tissueArea) + " um2"

// STEP 2: Load GeoJSON from Python
println "Step 2: Loading Python/LiverQuant GeoJSON..."
def geoJsonFile = new File(GEOJSON_FOLDER,
    slideName.replace('.svs', '_detections.geojson'))

def fatObjects = []
if (!geoJsonFile.exists()) {
    println "  WARNING: GeoJSON not found: " + geoJsonFile.getName()
    println "  Falling back to direct vacuole detection..."
    fatObjects = null  // signal to detect directly
} else {
    try {
        def imported = PathIO.readObjects(geoJsonFile)
        // Keep only fat/perinuclear/macrovesicular objects — exclude tissue outlines
        fatObjects = imported.findAll {
            def cls = it.getPathClass()?.toString()?.toLowerCase() ?: ''
            def nm  = it.getName()?.toLowerCase() ?: ''
            cls.contains('fat') || cls.contains('perinuclear') ||
            cls.contains('macro') || nm.contains('fat') ||
            nm.contains('perinuclear') || nm.contains('macro')
        }
        println "  Loaded ${fatObjects.size()} fat globules from GeoJSON"
        if (fatObjects.isEmpty()) {
            println "  WARNING: No fat objects found in GeoJSON — check labels"
            fatObjects = null
        }
    } catch (Exception e) {
        println "  WARNING: Could not load GeoJSON: ${e.message}"
        fatObjects = null
    }
}

// Fallback: detect vacuoles directly if GeoJSON failed
if (fatObjects == null) {
    println "  Detecting vacuoles directly..."
    selectAnnotations()
    runPlugin('qupath.imagej.detect.cells.WatershedCellDetection',
        '{"detectionImageBrightfield": "Optical density sum",' +
        '"requestedPixelSizeMicrons": 0.5,' +
        '"backgroundRadiusMicrons": 10.0,' +
        '"backgroundByReconstruction": true,' +
        '"medianRadiusMicrons": 0.0,' +
        '"sigmaMicrons": 2.0,' +
        '"minAreaMicrons": 20.0,' +
        '"maxAreaMicrons": 5000.0,' +
        '"threshold": 0.05,' +
        '"maxBackground": 0.5,' +
        '"watershedPostProcess": true,' +
        '"cellExpansionMicrons": 0.0,' +
        '"includeNuclei": false,' +
        '"smoothBoundaries": true,' +
        '"makeMeasurements": true}')
    fatObjects = getDetectionObjects().findAll { det ->
        def roi = det.getROI()
        if (roi == null) return false
        double a = roi.getArea()
        double p = roi.getLength()
        double c = (p > 0) ? (4 * Math.PI * a) / (p * p) : 0
        return c >= 0.5
    }
    println "  Vacuoles detected directly: ${fatObjects.size()}"
    clearDetections()
}

// STEP 3: Detect nuclei
println "Step 3: Detecting nuclei..."
selectAnnotations()
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection',
    '{"detectionImageBrightfield": "Hematoxylin OD",' +
    '"requestedPixelSizeMicrons": 0.5,' +
    '"backgroundRadiusMicrons": 8.0,' +
    '"backgroundByReconstruction": true,' +
    '"medianRadiusMicrons": 0.0,' +
    '"sigmaMicrons": 1.5,' +
    '"minAreaMicrons": 15.0,' +
    '"maxAreaMicrons": 400.0,' +
    '"threshold": 0.1,' +
    '"maxBackground": 2.0,' +
    '"watershedPostProcess": true,' +
    '"cellExpansionMicrons": 0.0,' +
    '"includeNuclei": false,' +
    '"smoothBoundaries": true,' +
    '"makeMeasurements": true}')

def nuclei = getDetectionObjects().collect()
int nucleiCount = nuclei.size()
println "  Nuclei detected: " + nucleiCount
if (nucleiCount == 0) {
    println "ERROR: No nuclei found. Skipping."
    return
}

// STEP 4: Build spatial grid index over fat globules
println "Step 4: Building spatial index over fat globules..."
double expansionPx = (!Double.isNaN(pixelSize) && pixelSize > 0) ?
    EXPANSION_MICRONS / pixelSize : EXPANSION_MICRONS
double gridSize = expansionPx * 5

def grid = [:]
fatObjects.each { fat ->
    def roi = fat.getROI()
    if (roi == null) return
    int gx = (int)(roi.getCentroidX() / gridSize)
    int gy = (int)(roi.getCentroidY() / gridSize)
    def key = "${gx}_${gy}"
    if (!grid.containsKey(key)) grid[key] = []
    grid[key] << fat
}
println "  Spatial grid: ${grid.size()} cells, ${fatObjects.size()} globules"

// STEP 5: Fast perinuclear matching
println "Step 5: Matching globules to nuclei..."
int totalDroplets    = 0
double totalLipidArea = 0
int nucleiWithLipid  = 0
int processed        = 0

nuclei.each { nucleus ->
    def nROI = nucleus.getROI()
    if (nROI == null) return

    double cx = nROI.getCentroidX()
    double cy = nROI.getCentroidY()
    def expandedGeom = nROI.getGeometry().buffer(expansionPx)

    int gxMin = (int)((cx - expansionPx) / gridSize) - 1
    int gxMax = (int)((cx + expansionPx) / gridSize) + 1
    int gyMin = (int)((cy - expansionPx) / gridSize) - 1
    int gyMax = (int)((cy + expansionPx) / gridSize) + 1

    int count   = 0
    double area = 0

    for (int gx = gxMin; gx <= gxMax; gx++) {
        for (int gy = gyMin; gy <= gyMax; gy++) {
            def nearby = grid["${gx}_${gy}"]
            if (nearby == null) continue
            nearby.each { fat ->
                def fROI = fat.getROI()
                if (fROI == null) return
                if (expandedGeom.intersects(fROI.getGeometry())) {
                    count++
                    double a = fROI.getArea()
                    area += (!Double.isNaN(pixelSize) && pixelSize > 0) ?
                        a * pixelSize * pixelSize : a
                }
            }
        }
    }

    totalDroplets   += count
    totalLipidArea  += area
    if (count > 0) nucleiWithLipid++

    processed++
    if (processed % 1000 == 0) {
        println "  Progress: ${processed}/${nucleiCount} nuclei done..."
    }
}

// STEP 6: Statistics
double pctWithLipid = (nucleiCount > 0) ? (nucleiWithLipid / nucleiCount) * 100 : 0
double avgDroplets  = (nucleiCount > 0) ? (double) totalDroplets / nucleiCount : 0
double avgArea      = (nucleiCount > 0) ? totalLipidArea / nucleiCount : 0
double avgSize      = (totalDroplets > 0) ? totalLipidArea / totalDroplets : 0

println "\n--- RESULTS ---"
println "Total nuclei:                  " + nucleiCount
println "Nuclei with perinuclear lipid: " + nucleiWithLipid + " (" + String.format("%.1f", pctWithLipid) + "%)"
println "Total perinuclear droplets:    " + totalDroplets
println "Avg droplets per nucleus:      " + String.format("%.2f", avgDroplets)
println "Avg lipid area per nucleus:    " + String.format("%.2f", avgArea) + " um2"
println "Avg individual droplet size:   " + String.format("%.2f", avgSize) + " um2"
println "Tissue area:                   " + String.format("%.1f", tissueArea) + " um2"
println "GeoJSON used:                  " + (geoJsonFile.exists() ? "Yes" : "No - direct detection")

// STEP 7: Save to CSV
def outFile = new File(outputCSV)
boolean needHeader = !outFile.exists() || outFile.length() == 0
outFile.withWriterAppend { w ->
    if (needHeader) {
        w.writeLine(
            "Image Name,Total Nuclei,Nuclei With Perinuclear Lipid," +
            "% Nuclei With Lipid,Total Perinuclear Droplets," +
            "Avg Droplets per Nucleus,Total Perinuclear Lipid Area (um2)," +
            "Avg Lipid Area per Nucleus (um2),Avg Droplet Size (um2)," +
            "Tissue Area (um2),GeoJSON Used"
        )
    }
    w.writeLine([
        slideName,
        nucleiCount,
        nucleiWithLipid,
        String.format("%.1f", pctWithLipid),
        totalDroplets,
        String.format("%.2f", avgDroplets),
        String.format("%.1f", totalLipidArea),
        String.format("%.2f", avgArea),
        String.format("%.2f", avgSize),
        String.format("%.1f", tissueArea),
        geoJsonFile.exists() ? "Yes" : "No"
    ].join(","))
}
println "Saved to: " + outputCSV
println "=== Done: " + slideName + " ==="
