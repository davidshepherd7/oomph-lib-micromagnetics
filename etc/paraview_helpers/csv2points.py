
from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()

TableToPoints1 = TableToPoints()
TableToPoints1.XColumn = 'x0'
TableToPoints1.YColumn = 'x1'
TableToPoints1.a2DPoints = 1

# SpreadSheetView1 = GetRenderView()
# DataRepresentation2 = Show()
# DataRepresentation2.FieldAssociation = 'Point Data'

# ltes = GetActiveSource()
# DataRepresentation1 = GetDisplayProperties(ltes)
# DataRepresentation1.Visibility = 0

# AnimationScene1 = GetAnimationScene()
# RenderView2 = CreateRenderView()
# RenderView2.CompressorConfig = 'vtkSquirtCompressor 0 3'
# RenderView2.UseLight = 1
# RenderView2.CameraPosition = [5.0, 5.0, 27.320508075688775]
# RenderView2.LightSwitch = 0
# RenderView2.OrientationAxesVisibility = 0
# RenderView2.CameraClippingRange = [27.047302994931886, 27.730315696824107]
# RenderView2.ViewTime = 0.0
# RenderView2.RemoteRenderThreshold = 3.0
# RenderView2.Background = [0.31999694819562063, 0.3400015259021897, 0.4299992370489052]
# RenderView2.CameraFocalPoint = [5.0, 5.0, 0.0]
# RenderView2.CameraParallelScale = 7.0710678118654755
# RenderView2.CenterOfRotation = [5.0, 5.0, 0.0]

# a1_error2_PVLookupTable = GetLookupTableForArray( " error2", 1, RGBPoints=[0.282843, 0.23, 0.299, 0.754, 0.282843, 0.706, 0.016, 0.15], VectorMode='Magnitude', NanColor=[0.25, 0.0, 0.0], ColorSpace='Diverging', ScalarRangeInitialized=1.0, AllowDuplicateScalars=1 )

# a1_error2_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )

# DataRepresentation3 = Show()
# DataRepresentation3.EdgeColor = [0.0, 0.0, 0.0]
# DataRepresentation3.PointSize = 5.0
# DataRepresentation3.SelectionPointFieldDataArrayName = ' '
# DataRepresentation3.ColorArrayName = ('POINT_DATA', ' error2')
# DataRepresentation3.LookupTable = a1_error2_PVLookupTable
# DataRepresentation3.Representation = 'Surface'

# AnimationScene1.ViewModules = [ SpreadSheetView1, RenderView2 ]

# a1_error2_PVLookupTable.ScalarOpacityFunction = a1_error2_PiecewiseFunction

# Render()
