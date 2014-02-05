try: paraview.simple
except: from paraview.simple import *

from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()

soln_pvd = GetActiveSource()
Calculator1 = Calculator()

RenderView1 = GetRenderView()
DataRepresentation1 = GetDisplayProperties(soln_pvd)
DataRepresentation2 = Show()
DataRepresentation2.ScaleFactor = 0.30000000000000004
DataRepresentation2.ScalarOpacityUnitDistance = 0.19840824707012364
DataRepresentation2.SelectionPointFieldDataArrayName = 'mvec'
DataRepresentation2.EdgeColor = [0.0, 0.0, 0.5000076295109483]

Calculator1.Function = 'V6*iHat+V7*jHat+V8*kHat'
Calculator1.ResultArrayName = 'mvec'

DataRepresentation1.Visibility = 0

Glyph1 = Glyph( GlyphType="Arrow", GlyphTransform="Transform2" )

Glyph1.Scalars = ['POINTS', 'V1']
Glyph1.SetScaleFactor = 0.30000000000000004
Glyph1.Vectors = ['POINTS', 'mvec']
Glyph1.GlyphTransform = "Transform2"
Glyph1.GlyphType = "Arrow"

Glyph1.SetScaleFactor = 0.29999999830887

DataRepresentation3 = Show()
DataRepresentation3.ScaleFactor = 0.32957027680240575
DataRepresentation3.SelectionPointFieldDataArrayName = 'V1'
DataRepresentation3.EdgeColor = [0.0, 0.0, 0.5000076295109483]

a1_V1_PVLookupTable = GetLookupTableForArray( "V1", 1, RGBPoints=[-9.850857734680176, 0.23, 0.299, 0.754, 5.319732189178467, 0.706, 0.016, 0.15], VectorMode='Magnitude', NanColor=[0.25, 0.0, 0.0], ColorSpace='Diverging', ScalarRangeInitialized=1.0, AllowDuplicateScalars=1 )

a1_V1_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )

DataRepresentation3.ColorArrayName = ('POINT_DATA', 'V1')
DataRepresentation3.LookupTable = a1_V1_PVLookupTable

a1_V1_PVLookupTable.ScalarOpacityFunction = a1_V1_PiecewiseFunction

a3_GlyphVector_PVLookupTable = GetLookupTableForArray( "GlyphVector", 3, RGBPoints=[1.0000000056370997, 0.23, 0.299, 0.754, 1.0000000056370997, 0.706, 0.016, 0.15], VectorMode='Magnitude', NanColor=[0.25, 0.0, 0.0], ColorSpace='Diverging', ScalarRangeInitialized=1.0, AllowDuplicateScalars=1 )

a3_GlyphVector_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )

RenderView1.CameraClippingRange = [4.5107056062273845, 8.820939956461126]

DataRepresentation3.ColorArrayName = ('POINT_DATA', 'GlyphVector')
DataRepresentation3.LookupTable = a3_GlyphVector_PVLookupTable

a3_GlyphVector_PVLookupTable.ScalarOpacityFunction = a3_GlyphVector_PiecewiseFunction
