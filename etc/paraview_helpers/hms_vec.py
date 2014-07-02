try: paraview.simple
except: from paraview.simple import *

from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()

soln_pvd = GetActiveSource()
Calculator2 = Calculator()

RenderView1 = GetRenderView()
DataRepresentation1 = GetDisplayProperties(soln_pvd)
DataRepresentation4 = Show()
DataRepresentation4.ScaleFactor = 0.30000000000000004
DataRepresentation4.ScalarOpacityUnitDistance = 0.19840824707012364
DataRepresentation4.SelectionPointFieldDataArrayName = 'Result'
DataRepresentation4.EdgeColor = [0.0, 0.0, 0.5000076295109483]

Calculator2.Function = 'V1*iHat+V2*jHat+V3*kHat'

DataRepresentation1.Visibility = 0

Glyph2 = Glyph( GlyphType="Arrow", GlyphTransform="Transform2" )

Glyph2.Scalars = ['POINTS', 'V1']
Glyph2.SetScaleFactor = 0.30000000000000004
Glyph2.Vectors = ['POINTS', 'Result']
Glyph2.GlyphTransform = "Transform2"
Glyph2.GlyphType = "Arrow"

Glyph2.SetScaleFactor = 0.0141834026997381

DataRepresentation5 = Show()
DataRepresentation5.ScaleFactor = 0.307979391515255
DataRepresentation5.SelectionPointFieldDataArrayName = 'V1'
DataRepresentation5.EdgeColor = [0.0, 0.0, 0.5000076295109483]

a1_V1_PVLookupTable = GetLookupTableForArray( "V1", 1, RGBPoints=[-16.631860733032227, 0.23, 0.299, 0.754, 6.445072174072266, 0.706, 0.016, 0.15] )

DataRepresentation5.ColorArrayName = ('POINT_DATA', 'V1')
DataRepresentation5.LookupTable = a1_V1_PVLookupTable

a3_GlyphVector_PVLookupTable = GetLookupTableForArray( "GlyphVector", 3, RGBPoints=[0.025796992328870942, 0.23, 0.299, 0.754, 21.1514829234553, 0.706, 0.016, 0.15] )

RenderView1.CameraClippingRange = [9.454775690340611, 14.731091870047191]

DataRepresentation5.ColorArrayName = ('POINT_DATA', 'GlyphVector')
DataRepresentation5.LookupTable = a3_GlyphVector_PVLookupTable

Render()
