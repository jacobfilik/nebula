package org.eclipse.nebular.visualisation.xygraph.plotlib;

import org.eclipse.draw2d.Graphics;
import org.eclipse.draw2d.PolygonDecoration;
import org.eclipse.draw2d.geometry.Point;
import org.eclipse.draw2d.geometry.PointList;
import org.eclipse.nebula.visualization.xygraph.figures.Axis;
import org.eclipse.swt.graphics.Color;

public class VectorField extends AbstractPlotArtist {

	private String         name;
	private String         dataName;
    private int[]          circleColor= new int[]{0,0,0};

	private int maximumArrowSize = 20;

	private double[] magnitude;
	private double[] orientation;
	private double[] x;
	private double[] y;
	private double maxMag = 1;
	private double minMag = 0;
//	private IDataset       vectors;
//	private IDataset       normalizedMagnitude;
//	private IDataset       normalizedAngle;
//	private List<IDataset> axes;
	private Axis xAxis;
	private Axis yAxis;
	private PolygonDecoration polyline;
//    private Map<RGB, Color> colorMap;
	private Color color;

	public VectorField() {

		this.polyline = new PolygonDecoration();

        PointList pl = new PointList();

        pl.addPoint(0, 0);
        pl.addPoint(-3, 1);
        pl.addPoint(-3, -1);
        polyline.setTemplate(pl);

		add(polyline);

//		colorMap = new HashMap<RGB, Color>(3);
	}

	public void dispose() {
//		// TODO Check this is called
//		for (Color col : colorMap.values()) {
//			col.dispose();
//		}
//		colorMap.clear();
	}

	/**
	 * This figure simply paints the vectors.
	 *
	 * TODO Downsampling!
	 */
	public void paint(Graphics graphics) {
		if (getLocalBackgroundColor() != null)
			graphics.setBackgroundColor(getLocalBackgroundColor());
		if (getLocalForegroundColor() != null)
			graphics.setForegroundColor(getLocalForegroundColor());

		graphics.pushState();
		try {
			paintArrows(graphics);
		} finally {
			graphics.popState();
		}
	}
	
	public void setColor(Color color) {
		this.color = color;
	}

	private void paintArrows(Graphics graphics) {

		for (int i = 0; i < magnitude.length; i++) {
			final double mag       = maximumArrowSize*(magnitude[i]-minMag)/(maxMag-minMag);
			final double angle     = orientation[i];
			final double xloc      = x[i];
			final double yloc      = y[i];

			// TODO Make normalize magnitude
			graphics.setForegroundColor(color);
			graphics.setBackgroundColor(color);
			paintArrow(graphics, xloc, yloc, mag, angle);
		}
	}

	/**
	 *
	 * @param xloc  - location of center of arrow
	 * @param yloc  - location of center of arrow
	 * @param length  - size of arrow, normalized from magnitude of vector.
	 * @param theta - anti-clockwise angle from 12 O'clock in radians.
	 */
	private void paintArrow(Graphics graphics, double xloc, double yloc, double length, double theta) {

		if (!xAxis.getRange().inRange(xloc) && !yAxis.getRange().inRange(yloc)) {
			return; // Nothing to draw
		}
		final int x    = xAxis.getValuePosition(xloc, false);
		final int y    = yAxis.getValuePosition(yloc, false);
		final double l = length/2d;

		
		
			final int xD = (int)Math.round(l*Math.sin(theta));
			final int yD = (int)Math.round(l*Math.cos(theta));
			final Point one = new Point(x-xD, y-yD);
			final Point two = new Point(x+xD, y+yD);
		


		graphics.drawLine(one, two);
		polyline.setLocation(one);
		polyline.setReferencePoint(two);
		polyline.setScale(l/3d, l/3d);

		polyline.paintFigure(graphics);
	}
	
	public void setData(double[] magnitude, double[] orientation, double[] x, double[] y) {
		this.magnitude = magnitude;
		this.orientation = orientation;
		this.x = x;
		this.y = y;
		
		maxMag = Double.NEGATIVE_INFINITY;
		minMag = Double.POSITIVE_INFINITY;
		
		for (int i = 0; i < magnitude.length; i++) {
			double val = magnitude[i];
			if (val > maxMag) maxMag = val;
			
			if (val < minMag) minMag = val;
		}
	}



	@Override
	protected void onAxesSet() {
		this.xAxis = axes.getPrimaryXAxis();
		this.yAxis = axes.getPrimaryYAxis();
		
	}

}
