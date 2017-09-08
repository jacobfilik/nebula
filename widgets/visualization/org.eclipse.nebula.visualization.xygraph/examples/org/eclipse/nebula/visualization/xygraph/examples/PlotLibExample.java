package org.eclipse.nebula.visualization.xygraph.examples;

import org.eclipse.draw2d.LightweightSystem;
import org.eclipse.nebula.visualization.xygraph.dataprovider.CircularBufferDataProvider;
import org.eclipse.nebula.visualization.xygraph.figures.DAxesFactory;
import org.eclipse.nebula.visualization.xygraph.figures.IXYGraph;
import org.eclipse.nebula.visualization.xygraph.figures.ToolbarArmedXYGraph;
import org.eclipse.nebula.visualization.xygraph.figures.Trace;
import org.eclipse.nebula.visualization.xygraph.figures.XYGraph;
import org.eclipse.nebular.visualisation.xygraph.plotlib.Axes;
import org.eclipse.nebular.visualisation.xygraph.plotlib.AxesImage;
import org.eclipse.nebular.visualisation.xygraph.plotlib.Line2D;
import org.eclipse.nebular.visualisation.xygraph.plotlib.VectorField;
import org.eclipse.nebula.visualization.xygraph.figures.Trace.PointStyle;
import org.eclipse.swt.widgets.Display;
import org.eclipse.swt.widgets.Shell;

public class PlotLibExample {
	public static void main(String[] args) {
		final Shell shell = new Shell();
		shell.setSize(300, 250);
		shell.open();

		// use LightweightSystem to create the bridge between SWT and draw2D
		final LightweightSystem lws = new LightweightSystem(shell);
		
		// create a new XY Graph.
		Axes axes = new Axes(new DAxesFactory());
		axes.setTitle("Simple Example");
		// set it as the content of LightwightSystem
		ToolbarArmedXYGraph toolbarArmedXYGraph = new ToolbarArmedXYGraph(axes);
		lws.setContents(toolbarArmedXYGraph);
		// create a trace data provider, which will provide the data to the
		// trace.
		CircularBufferDataProvider traceDataProvider = new CircularBufferDataProvider(false);
		traceDataProvider.setBufferSize(100);
		traceDataProvider.setCurrentXDataArray(new double[] { 10, 23, 34, 45, 56, 78, 88, 99 });
		traceDataProvider.setCurrentYDataArray(new double[] { 11, 44, 55, 45, 88, 98, 52, 23 });

		// create the trace
		Line2D line = new Line2D();
		line.setDataProvider(traceDataProvider);
		
//		Line2D trace = new Line2D("Trace1-XY Plot", xyGraph.getPrimaryXAxis(), xyGraph.getPrimaryYAxis(),
//				traceDataProvider);

		// set trace property
		line.setPointStyle(PointStyle.XCROSS);

		// add the trace to xyGraph
		axes.addArtist(line);
		
		VectorField f = new VectorField();
		
		f.setData(new double[] { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 },
				new double[] { 10, 20, 30, 40, 50, 60, 70, 80 }, 
				new double[] { 5, 15, 20, 25, 30, 35, 40, 45 }, 
				new double[] { 10, 10, 10, 10, 10, 10, 10, 10 });

		axes.addArtist(f);
		
		int[] shape = new int[] {10,10};
		
		double[] im = new double[shape[0]*shape[1]];
		for (int i = 0; i < im.length; i++) {
			im[i] = i;
		}
		
		
		
		AxesImage image = new AxesImage();
		image.setData(im, shape, new double[]{10,50,10,50});
		axes.addArtist(image);
		
//		axes.getPrimaryYAxis().setInverted(true);
		
		Display display = Display.getDefault();
		while (!shell.isDisposed()) {
			if (!display.readAndDispatch())
				display.sleep();
		}

	}
}
