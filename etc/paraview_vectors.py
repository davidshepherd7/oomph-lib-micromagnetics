try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

soln = GetActiveSource()
Calculator2 = Calculator()

Calculator2.AttributeMode = 'point_data'

Calculator2.Function = 'V1*iHat + V2*jHat + V3*kHat'

DataRepresentation2 = Show()

RenderView1 = GetRenderView()
DataRepresentation1 = GetDisplayProperties(soln)
a1_V4_PVLookupTable = GetLookupTableForArray( "V4", 1 )

DataRepresentation2.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation2.ColorAttributeType = 'POINT_DATA'
DataRepresentation2.ScalarOpacityFunction = []
DataRepresentation2.ColorArrayName = 'V4'
DataRepresentation2.ScalarOpacityUnitDistance = 0.2782305750496691
DataRepresentation2.Texture = []
DataRepresentation2.LookupTable = a1_V4_PVLookupTable

DataRepresentation1.Visibility = 0

Glyph1 = Glyph( GlyphType="Arrow" )

Glyph1.Scalars = ['POINTS', 'V1']
Glyph1.SetScaleFactor = 0.2
Glyph1.Vectors = ['POINTS', 'Result']
Glyph1.GlyphType = "Arrow"

Glyph1.SetScaleFactor = 0.395255568588082

DataRepresentation3 = Show()

a1_V1_PVLookupTable = GetLookupTableForArray( "V1", 1 )

DataRepresentation3.EdgeColor = [0.0, 0.0, 0.5000076295109483]

a1_V1_PiecewiseFunction = CreatePiecewiseFunction()

a3_GlyphVector_PVLookupTable = GetLookupTableForArray( "GlyphVector", 3, RGBPoints=[-1.0, 0.23, 0.299, 0.754, 1.0, 0.706, 0.016, 0.15], VectorMode='Magnitude', ColorSpace='Diverging', LockScalarRange=1 )

a3_GlyphVector_PiecewiseFunction = CreatePiecewiseFunction()

DataRepresentation2.Visibility = 0

Glyph1.MaximumNumberofPoints = 500

DataRepresentation3.ColorArrayName = 'GlyphVector'
DataRepresentation3.LookupTable = a3_GlyphVector_PVLookupTable
DataRepresentation3.Texture = []
DataRepresentation3.ColorAttributeType = 'POINT_DATA'

a1_V1_PVLookupTable.RGBPoints = [-1.0, 0.23, 0.299, 0.754, 1.0, 0.706, 0.016, 0.15]
a1_V1_PVLookupTable.VectorMode = 'Magnitude'
a1_V1_PVLookupTable.ColorSpace = 'Diverging'
a1_V1_PVLookupTable.LockScalarRange = 1

RenderView1.CameraViewUp = [-0.04477112623956685, -0.9917307302341752, 0.12027345910229019]
RenderView1.CameraPosition = [4.307604873906881, -0.009073836104138108, 1.5286624402970659]
RenderView1.CameraClippingRange = [1.899628019870475, 8.15553918000676]

Glyph1.MaximumNumberofPoints = 1000

RenderView1.CameraViewUp = [-0.25137000935902554, -0.9582983952798982, 0.1359312473231053]
RenderView1.CameraPosition = [1.916237570151724, 0.08585492340087403, 4.148856884426094]
RenderView1.CameraClippingRange = [1.784159586143178, 8.18493042329335]

Render()
