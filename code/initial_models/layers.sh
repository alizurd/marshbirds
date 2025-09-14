—check and see if there are 4 bands
gdalinfo ~/Desktop/naip/tiffs/bx_clip.tif

—extract each band and export as a new layer to desktop 
gdal_translate -b 3 /Users/alyssabueno/Desktop/m_4007310_ne_18_060_20211105.tif ~/Desktop/three_band.tif