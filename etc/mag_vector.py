try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

DataRepresentation1 = Show()

soln = GetActiveSource()
RenderView1 = GetRenderView()
DataRepresentation1.ScalarOpacityUnitDistance = 0.8039467687661698
DataRepresentation1.EdgeColor = [0.0, 0.0, 0.5000076295109483]

Calculator1 = Calculator()

Calculator1.AttributeMode = 'point_data'

Calculator1.Function = 'V6*iHat+V7*jHat+V8*kHat'

DataRepresentation2 = Show()

DataRepresentation1.Visibility = 0

DataRepresentation2.ScalarOpacityUnitDistance = 0.8039467687661698
DataRepresentation2.Texture = []
DataRepresentation2.EdgeColor = [0.0, 0.0, 0.5000076295109483]

Glyph1 = Glyph( GlyphType="Arrow" )

Glyph1.Scalars = ['POINTS', 'V1']
Glyph1.SetScaleFactor = 0.2
Glyph1.Vectors = ['POINTS', 'Result']
Glyph1.GlyphType = "Arrow"

DataRepresentation3 = Show()

a1_V1_PVLookupTable = GetLookupTableForArray( "V1", 1 )

DataRepresentation3.EdgeColor = [0.0, 0.0, 0.5000076295109483]

a1_V1_PiecewiseFunction = CreatePiecewiseFunction()

a3_GlyphVector_PVLookupTable = GetLookupTableForArray( "GlyphVector", 3, RGBPoints=[1.0, 0.23, 0.299, 0.754, 1.0, 0.706, 0.016, 0.15], VectorMode='Magnitude', ColorSpace='Diverging', ScalarRangeInitialized=1.0, LockScalarRange=1 )

a3_GlyphVector_PiecewiseFunction = CreatePiecewiseFunction()

ScalarBarWidgetRepresentation1 = CreateScalarBar( Title='GlyphVector Magnitude', Position2=[0.13, 0.5], Enabled=1, LabelFontSize=12, LookupTable=a3_GlyphVector_PVLookupTable, TitleFontSize=12, Position=[0.87, 0.25] )
GetRenderView().Representations.append(ScalarBarWidgetRepresentation1)

RenderView1.CameraViewUp = [-0.002002264440739648, 0.9990724496133485, -0.04301431576455075]
RenderView1.CameraPosition = [1.3545918110218906, 0.28460896512998785, 6.547418457171226]
RenderView1.CameraClippingRange = [4.135144960796614, 9.972760497248697]

DataRepresentation2.Visibility = 0

DataRepresentation3.ColorArrayName = 'GlyphVector'
DataRepresentation3.LookupTable = a3_GlyphVector_PVLookupTable
DataRepresentation3.Texture = []
DataRepresentation3.ColorAttributeType = 'POINT_DATA'

a1_V1_PVLookupTable.RGBPoints = [-1.0, 0.23, 0.299, 0.754, 1.0, 0.706, 0.016, 0.15]
a1_V1_PVLookupTable.VectorMode = 'Magnitude'
a1_V1_PVLookupTable.ColorSpace = 'Diverging'
a1_V1_PVLookupTable.LockScalarRange = 1

Render()
