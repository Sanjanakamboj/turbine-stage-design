# ParaView batch script — CFD colormaps from flow.vtu
# Run:  /Applications/ParaView-6.1.0.app/Contents/bin/pvbatch paraview_colormaps.py
from paraview.simple import *
import os

paraview.simple._DisableFirstRenderCameraReset()
W = '/Users/sanju/Projects/Python/Turbine Stage Design'
R = os.path.join(W, 'Results')

flow = XMLUnstructuredGridReader(FileName=[os.path.join(W, 'flow.vtu')])
flow.UpdatePipeline()

# entropy rise:  ds = cp*ln(T/Tt0) - Rgas*ln(P/Pt0)
entropy = Calculator(Input=flow)
entropy.ResultArrayName = 'EntropyRise'
entropy.AttributeType = 'Point Data'
entropy.Function = '1156.45*ln(Temperature/1578.0609)-287.0*ln(Pressure/1521338.145)'
entropy.UpdatePipeline()

view = CreateView('RenderView')
view.ViewSize = [850, 1400]
view.OrientationAxesVisibility = 0
view.CameraParallelProjection = 1
view.UseColorPaletteForBackground = 0
view.BackgroundColorMode = 'Single Color'
view.Background = [1.0, 1.0, 1.0]

xmin, xmax, ymin, ymax = flow.GetDataInformation().GetBounds()[:4]
xc, yc = 0.5 * (xmin + xmax), 0.5 * (ymin + ymax)


def render_field(source, field, comp, title, preset, fname):
    disp = Show(source, view)
    disp.Representation = 'Surface'
    color = ('POINTS', field) if comp is None else ('POINTS', field, comp)
    ColorBy(disp, color)
    lut = GetColorTransferFunction(field)
    try:
        lut.ApplyPreset(preset, True)
    except Exception as e:
        print('preset fallback for', field, ':', e)
    disp.RescaleTransferFunctionToDataRange(True, False)
    disp.SetScalarBarVisibility(view, True)
    sb = GetScalarBar(lut, view)
    sb.Title = title
    sb.ComponentTitle = ''
    sb.TitleColor = [0, 0, 0]
    sb.LabelColor = [0, 0, 0]
    sb.TitleFontSize = 16
    sb.LabelFontSize = 14
    sb.AutoOrient = 0
    sb.Orientation = 'Vertical'
    sb.WindowLocation = 'Any Location'
    sb.Position = [0.80, 0.32]
    sb.ScalarBarLength = 0.42
    sb.ScalarBarThickness = 18
    # 2-D top-down view (look along -z)
    view.CameraPosition = [xc, yc, 1.0]
    view.CameraFocalPoint = [xc, yc, 0.0]
    view.CameraViewUp = [0, 1, 0]
    view.ResetCamera()
    Render()
    out = os.path.join(R, fname)
    SaveScreenshot(out, view, ImageResolution=[850, 1400], TransparentBackground=0)
    print('saved', out)
    Hide(source, view)


render_field(flow,    'Mach',        None,        'Mach number',            'Jet',                  'PV_Mach.png')
render_field(flow,    'Pressure',    None,        'Static Pressure [Pa]',   'Cool to Warm',         'PV_Pressure.png')
render_field(flow,    'Temperature', None,        'Temperature [K]',        'Inferno',              'PV_Temperature.png')
render_field(flow,    'Velocity',    'Magnitude', 'Velocity [m/s]',         'Viridis',              'PV_Velocity.png')
render_field(entropy, 'EntropyRise', None,        'Entropy rise [J/kg/K]',  'Black-Body Radiation', 'PV_Entropy.png')
print('DONE')
