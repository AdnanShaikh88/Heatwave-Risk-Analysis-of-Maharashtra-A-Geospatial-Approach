// ============================================================================
// HEATWAVE RISK MAPPING - SMOOTH 1KM VERSION
// Interpolates coarse hazard data to match fine exposure/vulnerability
// ============================================================================

var aoi = ee.FeatureCollection("FAO/GAUL/2015/level1")
            .filter(ee.Filter.eq('ADM1_NAME', 'Maharashtra'))
            .geometry();
Map.centerObject(aoi, 7);

var baselineStart = '1991-01-01';
var baselineEnd = '2020-12-31';
var studyStart = '2023-03-01';
var studyEnd = '2023-05-31';

// ============================================================================
// STEP 1: CALCULATE HEAT INDEX
// ============================================================================

var era5 = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR')
  .filterBounds(aoi)
  .select(['temperature_2m', 'dewpoint_temperature_2m']);

var kToC = function(image) {
  return image.subtract(273.15).copyProperties(image, ['system:time_start']);
};

var calculateRH = function(image) {
  var t = image.select('temperature_2m');
  var td = image.select('dewpoint_temperature_2m');
  var es = t.expression('6.11 * pow(10, (7.5 * T) / (237.3 + T))', {'T': t});
  var e = td.expression('6.11 * pow(10, (7.5 * T) / (237.3 + T))', {'T': td});
  var rh = e.divide(es).multiply(100).clamp(0, 100);
  return image.addBands(rh.rename('RH'));
};

var calculateHI = function(image) {
  var t = image.select('temperature_2m');
  var rh = image.select('RH');
  
  var hi = t.expression(
    'c1 + (c2*T) + (c3*RH) + (c4*T*RH) + (c5*T*T) + (c6*RH*RH) + (c7*T*T*RH) + (c8*T*RH*RH) + (c9*T*T*RH*RH)', {
    'T': t, 'RH': rh,
    'c1': -8.78469475556, 'c2': 1.61139411, 'c3': 2.33854883889,
    'c4': -0.14611605, 'c5': -0.012308094, 'c6': -0.0164248277778,
    'c7': 0.002211732, 'c8': 0.00072546, 'c9': -0.000003582
  });
  
  // FIXED: Corrected parentheses
  var simpleHI = t.expression('0.5 * (T + 16.92 + abs(T - 16.92) + 0.18 * RH)', {'T': t, 'RH': rh});
  var finalHI = hi.where(t.lt(26.7), simpleHI);
  
  return image.addBands(finalHI.rename('HI'));
};

var hiCollection = era5.map(kToC).map(calculateRH).map(calculateHI);
var dailyHI = hiCollection.select('HI');

// ============================================================================
// STEP 2: CALCULATE 90TH PERCENTILE THRESHOLDS
// ============================================================================

var baselineHI = dailyHI.filterDate(baselineStart, baselineEnd);

function percentileByDecade(start, end) {
  var subset = baselineHI.filterDate(start, end);
  return ee.ImageCollection.fromImages(
    ee.List.sequence(1, 12).map(function(m) {
      var monthImgs = subset.filter(ee.Filter.calendarRange(m, m, 'month'));
      return monthImgs.reduce(ee.Reducer.percentile([90]))
                      .rename('HI_threshold')
                      .set('month', m);
    })
  );
}

var baseline1 = percentileByDecade('1991-01-01', '2000-12-31');
var baseline2 = percentileByDecade('2001-01-01', '2010-12-31');
var baseline3 = percentileByDecade('2011-01-01', '2020-12-31');

var merged = baseline1.merge(baseline2).merge(baseline3);

var dailyThresholds = ee.ImageCollection.fromImages(
  ee.List.sequence(1, 12).map(function(m) {
    var monthly = merged.filter(ee.Filter.eq('month', m));
    return monthly.mean().set('month', m);
  })
);

// ============================================================================
// STEP 3: IDENTIFY HEATWAVE DAYS
// ============================================================================

var studyHI = dailyHI.filterDate(studyStart, studyEnd);

var findHeatwaveDays = function(image) {
  var m = image.date().get('month');  
  var threshold = ee.Image(
    dailyThresholds.filter(ee.Filter.eq('month', m)).first()
  );
  return image.gt(threshold)
              .rename('is_heatwave')
              .copyProperties(image, ['system:time_start']);
};

var heatwaveBinary = studyHI.map(findHeatwaveDays);

// Sum heatwave days - THIS IS COARSE (~11km actual resolution)
var hazardImage_coarse = heatwaveBinary.reduce(ee.Reducer.sum())
                                       .clip(aoi)
                                       .unmask(0);

print('Original Hazard Resolution:', hazardImage_coarse.projection().nominalScale());

// ============================================================================
// STEP 4: REPROJECT HAZARD TO 1KM (SMOOTH INTERPOLATION)
// ============================================================================
// Change this line in STEP 4:
// Force hazard to 1km using bilinear interpolation
var hazardImage = hazardImage_coarse
  .resample('bicubic')  // Smooth interpolation
  .reproject({
    crs: 'EPSG:4326',
    scale: 1000  // Force to 1km
  });

print('New Hazard Resolution:', hazardImage.projection().nominalScale());

var hazardVis = {
  min: 0, max: 30, 
  palette: ['#ffffcc', '#ffeda0', '#fed976', '#feb24c', '#fd8d3c', '#fc4e2a', '#e31a1c', '#b10026']
};
Map.addLayer(hazardImage, hazardVis, 'Hazard: Heatwave Days (1km smooth)');

// ============================================================================
// STEP 5: EXPOSURE - POPULATION (ALREADY 1KM, NO CHANGES)
// ============================================================================

var gpw = ee.ImageCollection("CIESIN/GPWv411/GPW_Population_Density")
            .filterDate('2020-01-01', '2020-12-31')
            .first();
var exposureImage = gpw.clip(aoi).unmask(0); 

var popVis = {min: 0, max: 1000, palette: ['#ffffe5', '#fec44f', '#d95f0e']};
Map.addLayer(exposureImage, popVis, 'Exposure: Population (1km)', false);

print('Exposure Resolution:', exposureImage.projection().nominalScale());

// ============================================================================
// STEP 6: VULNERABILITY - COMPOSITE INDEX (ALREADY 1KM)
// ============================================================================

var ndvi = ee.ImageCollection('MODIS/061/MOD13A3')
            .filterDate(studyStart, studyEnd)
            .select('NDVI')
            .mean()
            .clip(aoi);
var v_env = ee.Image(1).subtract(ndvi.unitScale(0, 10000)).clamp(0, 1);

var nightLST = ee.ImageCollection('MODIS/061/MYD11A1')
            .filterDate(studyStart, studyEnd)
            .select('LST_Night_1km')
            .mean()
            .clip(aoi);
var nightLST_C = nightLST.multiply(0.02).subtract(273.15);

var nightStats = nightLST_C.reduceRegion({
  reducer: ee.Reducer.minMax(), 
  geometry: aoi, 
  scale: 1000, 
  bestEffort: true
});

var v_uhi = nightLST_C.unitScale(
  nightStats.get('LST_Night_1km_min'), 
  nightStats.get('LST_Night_1km_max')
).clamp(0, 1);

var vulnerabilityImage = v_env.multiply(0.4).add(v_uhi.multiply(0.6)).unmask(0);

var vulnVis = {min: 0, max: 1, palette: ['#f7f7f7', '#d95f0e']};
Map.addLayer(vulnerabilityImage, vulnVis, 'Vulnerability Index', false);

print('Vulnerability Resolution:', vulnerabilityImage.projection().nominalScale());

// ============================================================================
// STEP 7: NO HARMONIZATION NEEDED - ALL AT 1KM NOW!
// ============================================================================

// Just ensure they're all masked consistently
var exposureHarmonized = exposureImage.unmask(0);
var vulnerabilityHarmonized = vulnerabilityImage.unmask(0);

print('✓ All layers now at 1km - no reduceResolution needed!');

// ============================================================================
// STEP 8: NORMALIZE ALL LAYERS (0-1)
// ============================================================================

var getStats = function(image) {
  return image.reduceRegion({
    reducer: ee.Reducer.minMax(),
    geometry: aoi,
    scale: 1000,  // Use 1km for all stats
    maxPixels: 1e9,
    bestEffort: true
  });
};

var hazStats = getStats(hazardImage.rename('val'));
var expStats = getStats(exposureHarmonized.rename('val'));
var vulnStats = getStats(vulnerabilityHarmonized.rename('val'));

print('Hazard range:', hazStats);
print('Exposure range:', expStats);
print('Vulnerability range:', vulnStats);

var normalizedHazard = hazardImage.unitScale(
  ee.Number(hazStats.get('val_min')), 
  ee.Number(hazStats.get('val_max'))
).clamp(0, 1);

var normalizedExposure = exposureHarmonized.unitScale(
  ee.Number(expStats.get('val_min')), 
  ee.Number(expStats.get('val_max'))
).clamp(0, 1);

var normalizedVulnerability = vulnerabilityHarmonized.unitScale(
  ee.Number(vulnStats.get('val_min')), 
  ee.Number(vulnStats.get('val_max'))
).clamp(0, 1);


// ============================================================================
// STEP 9: CALCULATE FINAL RISK MAP (H × E × V)
// ============================================================================

// --- FIX: Rename the final band to 'risk' for consistent stats and exports
var finalRiskMap = normalizedHazard
  .multiply(normalizedExposure)
  .multiply(normalizedVulnerability)
  .clip(aoi)
  .rename('risk'); // <-- RENAMED HERE

// Check final stats
var riskStats = finalRiskMap.reduceRegion({ // Now runs on the 'risk' band
  reducer: ee.Reducer.minMax(),
  geometry: aoi,
  scale: 1000,
  maxPixels: 1e9,
  bestEffort: true
});
// --- END FIX

print('Final Risk range:', riskStats);

// Adjust visualization based on actual data range
// var riskMax = ee.Number(riskStats.get('constant_max')).getInfo(); // Old
var riskMax = ee.Number(riskStats.get('risk_max')).getInfo(); // <-- FIX: Get 'risk_max'
print('Suggested max for visualization:', riskMax);

var riskVis = {
  min: 0, 
  max: 0.1, // Your hard-coded value for good color contrast
  
  palette: ['#08306b', '#08519c', '#3182bd', '#fd8d3c', '#f03b20', '#bd0026', '#800026']


};
Map.addLayer(finalRiskMap, riskVis, 'Final Risk Map (H × E × V) - 1km');

// ============================================================================
// STEP 10: ADD LEGENDS
// ============================================================================

function createLegend(title, palette, min, max, position) {
  var legend = ui.Panel({
    style: {
      position: position,
      padding: '8px 15px',
      backgroundColor: 'rgba(255, 255, 255, 0.9)'
    }
  });
  
  var legendTitle = ui.Label({
    value: title,
    style: {fontWeight: 'bold', fontSize: '13px', margin: '0 0 4px 0'}
  });
  legend.add(legendTitle);

  var colorBar = ui.Thumbnail({
    image: ee.Image.pixelLonLat().select(0),
    params: {
      bbox: [0, 0, 1, 0.1],
      dimensions: '100x10',
      format: 'png',
      min: 0, max: 1,
      palette: palette,
    },
    style: {stretch: 'horizontal', margin: '0px 8px', maxHeight: '20px'},
  });
  legend.add(colorBar);
  
  var labels = ui.Panel({
    widgets: [
      ui.Label(min, {margin: '4px 8px', fontSize: '11px'}),
      ui.Label(((min + max) / 2).toFixed(1), 
               {margin: '4px 8px', textAlign: 'center', stretch: 'horizontal', fontSize: '11px'}),
      ui.Label(max, {margin: '4px 8px', fontSize: '11px'})
    ],
    layout: ui.Panel.Layout.flow('horizontal')
  });
  legend.add(labels);
  
  return legend;
}

Map.add(createLegend('Hazard: Days', hazardVis.palette, hazardVis.min, hazardVis.max, 'bottom-center'));
Map.add(createLegend('Exposure: Pop/km²', popVis.palette, popVis.min, popVis.max, 'bottom-left'));
Map.add(createLegend('Vulnerability', vulnVis.palette, vulnVis.min, vulnVis.max, 'top-right'));
Map.add(createLegend('Risk (H×E×V)', riskVis.palette, riskVis.min, riskVis.max, 'bottom-right'));

// ============================================================================
// STEP 11: DISTRICT-LEVEL ANALYSIS
// ============================================================================

var districts = ee.FeatureCollection("FAO/GAUL/2015/level2")
  .filter(ee.Filter.eq('ADM1_NAME', 'Maharashtra'));

var districtOutline = ee.Image().byte().paint({
  featureCollection: districts,
  color: 1,
  width: 1.5
});
Map.addLayer(districtOutline, {palette: '000000'}, 'District Boundaries');

// Calculate district statistics
var districtRisk = finalRiskMap.reduceRegions({
  collection: districts,
  reducer: ee.Reducer.mean().combine({
    reducer2: ee.Reducer.stdDev(),
    sharedInputs: true
  }).combine({
    reducer2: ee.Reducer.max(),
    sharedInputs: true
  }),
  scale: 1000,
  crs: 'EPSG:4326'
});

var riskTable = districtRisk.map(function(feature) {
  var name = feature.get('ADM2_NAME');
  var riskMean = feature.get('mean');
  var riskStd = feature.get('stdDev');
  var riskMax = feature.get('max');
  return ee.Feature(null, {
    'District': name, 
    'Risk_Mean': riskMean,
    'Risk_StdDev': riskStd,
    'Risk_Max': riskMax
  });
});

print('═══════════════════════════════════════════════════════════');
print('DISTRICT RISK RANKINGS (Descending by Mean Risk):');
print('═══════════════════════════════════════════════════════════');
print(riskTable.sort('Risk_Mean', false));

// ============================================================================
// STEP 12: EXPORT OUTPUTS
// ============================================================================

Export.table.toDrive({
  collection: riskTable.sort('Risk_Mean', false),
  description: 'Maharashtra_Heatwave_Risk_Districts_1km',
  fileFormat: 'CSV'
});

Export.image.toDrive({
  image: finalRiskMap,
  description: 'Maharashtra_Heatwave_Risk_Raw_1km',
  scale: 1000,
  region: aoi,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

Export.image.toDrive({
  image: finalRiskMap.visualize(riskVis),
  description: 'Maharashtra_Heatwave_Risk_Visualized_1km',
  scale: 1000,
  region: aoi,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

// ============================================================================
// COMPARISON: SHOW BEFORE/AFTER
// ============================================================================

Map.addLayer(hazardImage_coarse, hazardVis, 'Hazard (Original ~11km blocky)', false);
Map.addLayer(hazardImage, hazardVis, 'Hazard (Smoothed to 1km)', false);

print('');
print('✅ CODE EXECUTION COMPLETE');
print('');
print('📊 WHAT CHANGED:');
print('   • Hazard layer interpolated from ~11km → 1km');
print('   • All three layers now aligned at 1km resolution');
print('   • No reduceResolution() needed - simpler workflow');
print('   • Map should look complete and smooth now');
print('');
print('⚠️  IMPORTANT NOTE FOR YOUR REPORT:');
print('   "Hazard data (ERA5-Land) was bilinearly interpolated from');
print('   native ~11km resolution to 1km to match exposure and');
print('   vulnerability layers. This improves visualization but does');
print('   not add new information to the climate data."');
print('');
print('🔍 Check the print outputs above for actual data ranges');
print('📁 Check Tasks tab to run exports');
