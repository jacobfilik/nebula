package org.eclipse.nebular.visualisation.xygraph.plotlib;

import org.eclipse.draw2d.Graphics;
import org.eclipse.draw2d.PolygonDecoration;
import org.eclipse.draw2d.geometry.Point;
import org.eclipse.draw2d.geometry.PointList;
import org.eclipse.nebula.visualization.xygraph.figures.Axis;
import org.eclipse.nebula.visualization.xygraph.linearscale.Range;
import org.eclipse.swt.graphics.Color;
import org.eclipse.swt.graphics.Image;
import org.eclipse.swt.graphics.ImageData;
import org.eclipse.swt.graphics.PaletteData;
import org.eclipse.swt.widgets.Display;

public class AxesImage extends AbstractPlotArtist {

	private double[] magnitude;
	private double[] x;
	private double[] y;
	private double maxMag = 1;
	private double minMag = 0;
	private int currentDownSampleBin = 1;
//	private IDataset       vectors;
//	private IDataset       normalizedMagnitude;
//	private IDataset       normalizedAngle;
//	private List<IDataset> axes;
	private Axis xAxis;
	private Axis yAxis;
	private PolygonDecoration polyline;
	private PaletteData palette = new PaletteData(0xff, 0xff00, 0xff0000);
//    private Map<RGB, Color> colorMap;
	private Color color;

	public AxesImage() {

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
//		this.orientation = orientation;
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
	
private boolean buildImageRelativeToAxes(ImageData imageData) {
		
		// slice data to get current zoom area
		/**     
		 *      x1,y1--------------x2,y2
		 *        |                  |
		 *        |                  |
		 *        |                  |
		 *      x3,y3--------------x4,y4
		 */
		
		ImageData data = imageData;
		
		
		if (xAxis == null || yAxis == null) return false;
		
		
		boolean xDataInc = x[0] < x[x.length-1];
		boolean yDataInc = y[0] < y[y.length-1];
		
		
		//Get the axes coodinates visible on screen
		double[] da = getImageCoords(1, false, x, y);
		double minX = da[0];
		double minY = da[1];
		double maxX = da[2];
		double maxY = da[3];
		
		//Get the data points per axis step
		double xAxValPerPoint =  Math.abs(getAxisValuePerDataPoint(minX,maxX,xAxis));
		double yAxValPerPoint =  Math.abs(getAxisValuePerDataPoint(minY,maxY,yAxis));
		
		if (Double.isNaN(xAxValPerPoint)|| Double.isNaN(yAxValPerPoint)) {
			return false;
		}
		
		//get x and y start position in data array (floored)
		int xPix = getPositionInAxis(xDataInc ? minX : maxX,xAxis,true)/currentDownSampleBin;
		int yPix = getPositionInAxis(yDataInc ? minY : maxY,yAxis,true)/currentDownSampleBin;
		
		//determine sub pixel offset in data frame
		double xdiff = (xAxis.getDouble(xPix)-(xDataInc ? minX : maxX));
		double ydiff = (yAxis.getDouble(yPix)-(yDataInc ? minY : maxY));
		
		//Determine the corresponding number of data points in x and y (floor min, ceil max)
		int xDataPoints = getNumberOfDataPoints(minX,maxX,getAxes().get(0), xDataInc);
		int yDataPoints = getNumberOfDataPoints(minY,maxY,getAxes().get(1), yDataInc);
		
		int xDataPointsDS = xDataPoints/currentDownSampleBin;
		int yDataPointsDS = yDataPoints/currentDownSampleBin;
		
		//Get matching screen co-ordinates in pixels
		double[] screenCoords = getImageCoords(1, true, xAxis, yAxis);
		
		//calculate full number of data points
		double realx = (maxX-minX)/xAxValPerPoint;
		double realy = (maxY-minY)/yAxValPerPoint;
		
		// Ratio of screen pixels to full number of data points
		double xScale = Math.abs((screenCoords[2]-screenCoords[0]) / (realx));
		double yScale = Math.abs((screenCoords[3]-screenCoords[1]) / (realy));
		
		//FIXME Handle minimum cases
		
		// Size of image data in index
		final int xSize = imageData.width;
		final int ySize = imageData.height;
		
		double[] minMaxX = getMinMax(x);
		double[] minMaxY = getMinMax(y);
		
		double xmin = minMaxX[0];
		double xmax = minMaxX[1];
		double ymin = minMaxY[0];
		double ymax = minMaxY[1];
		
		double xLen = x.length;
		double yLen = y.length;
		
		double xp = ((xmax-xmin)/(xLen-1))/2;
		double yp = ((ymax-ymin)/(yLen-1))/2;
		
		double xOffset = (((-xdiff+xp)/xAxValPerPoint))*xScale;
		double yOffset = (((-ydiff+yp)/yAxValPerPoint))*yScale;
//		xOffset = 0;
		

		//FIXME origins
		
		int width  = xPix+xDataPointsDS > xSize  ? Math.min(xSize, xDataPointsDS) : xDataPointsDS;
		int height = yPix+yDataPointsDS > ySize  ? Math.min(ySize, yDataPointsDS) : yDataPointsDS;
		
		int scaleWidth  = Math.max(1, (int) (xDataPoints*xScale));
		int scaleHeight = Math.max(1, (int) (yDataPoints*yScale));
		
//		 Force a minimum size on the system
//		if (width <= MINIMUM_ZOOM_SIZE) {
//			if (width > imageData.width) width = MINIMUM_ZOOM_SIZE;
//			isMaximumZoom = true;
//		}
//		if (height <= MINIMUM_ZOOM_SIZE) {
//			if (height > imageData.height) height = MINIMUM_ZOOM_SIZE;
//			isMaximumZoom = true;
//		}
		
		try {
			// Slice the data.
			// Pixel slice on downsampled data = fast!
			data = sliceImageData(data, width, height, xPix, yPix, ySize);

			// create the scaled image
			// We are suspicious if the algorithm wants to create an image
			// bigger than the screen size and in that case do not scale
			// Fix to http://jira.diamond.ac.uk/browse/SCI-926
			boolean proceedWithScale = true;
			try {
				if (screenRectangle == null) {
					screenRectangle = Display.getCurrent().getPrimaryMonitor().getClientArea();
				}
				if (scaleWidth>screenRectangle.width*2 ||  scaleHeight>screenRectangle.height*2) {
					logger.error("Image scaling algorithm has malfunctioned and asked for an image bigger than the screen!");
					logger.debug("scaleWidth="+scaleWidth);
					logger.debug("scaleHeight="+scaleHeight);
					proceedWithScale = false;
				}
			} catch (Throwable ne) {
				proceedWithScale = true;
			}

			if (!xDataInc) scaleWidth*=-1;
			if (!yDataInc) scaleHeight*=-1;
			
			
			Image scaledImage = null;
			if (proceedWithScale) {
				data = data!=null ? data.scaledTo(scaleWidth, scaleHeight) : null;
				scaledImage = data!=null ? new Image(Display.getDefault(), data) : null;
			} else if (scaledImage==null) {
				scaledImage = data!=null ? new Image(Display.getDefault(), data) : null;
			}
			
			scaledData.setX(screenCoords[0]);
			scaledData.setY(screenCoords[1]);
			scaledData.setXoffset(xOffset);
			scaledData.setYoffset(yOffset);
			scaledData.setScaledImage(scaledImage);
			
		} catch (IllegalArgumentException ne) {
			
			logger.error("Image scaling has malfunctioned");
			logger.debug("Trace name = "+getName());
			logger.debug("scaleWidth = "+scaleWidth);
			logger.debug("scaleHeight = "+scaleHeight);
			logger.debug("width = "+width);
			logger.debug("height = "+height);
			logger.debug("xPix = "+xPix);
			logger.debug("yPix = "+yPix);
			
			throw ne;
		}
		
		
		return true;
		
	}

private final double getAxisValuePerDataPoint(double min, double max, double[] axes) {
	int minp = getPositionInAxis(min,axes);
	int maxp = getPositionInAxis(max,axes);
	
	return (axes[maxp]-axes[minp])/(maxp-minp);
	
}

private final int getPositionInAxis(double val, double[] axis) {
	if (axis.length == 1) return 0;
	
	int pos = 0;
	double v = axis[pos];
	if (val - v == 0) return 0;
	
	double minVal = Math.abs(val-v);
	
	for (int i = 1; i < axis.length; i++) {
		v = axis[i];
		
		double diff = v-val;
		
		if (diff == 0) return i;
		
		double temp =  Math.abs(val-v);
		
		if (temp < minVal) {
			minVal = temp;
			pos = i;
		}
		
	}
	
	return pos;
}

//private final int getPositionInAxis(double val, IDataset axis, boolean floor) {
//	if (axis.getSize() == 1) return 0;
//	int pos = 0;
//	double dif = 0;
//	double test = Double.MAX_VALUE;
//	Dataset a = DatasetUtils.convertToDataset(axis.clone());
//	IndexIterator it = a.getIterator();
//	while (it.hasNext()) {
//		double d = a.getDouble(it.index);
//		double ad = Math.abs(d-val);
//		if (ad < test) {
//			test = ad;
//			pos = it.index;
//			dif = d-val;
//		}
//	}
//	
//	if (floor) {
//		if (dif > 0) pos--;
//	} else {
//		if (dif < 0) pos++;
//	}
//	
//	return pos;
//	
//}
private ImageData sliceImageData(ImageData imageData, int width, int height, int xPix, int yPix, int ySize) {
	ImageData data;
	if (imageData.depth <= 8) {
		// NOTE Assumes 8-bit images
		final int size   = width*height;
		final byte[] pixels = new byte[size];
		for (int y = 0; y < height && (yPix+y)<ySize ; y++) {
			imageData.getPixels(xPix, yPix+y, width, pixels, width*y);
		}
		data = new ImageData(width, height, imageData.depth, palette, 1, pixels);
		if (imageData.alphaData != null) {
			final byte[] alphas = new byte[size];
			for (int y = 0; y < height && (yPix+y)<ySize ; y++) {
				imageData.getAlphas(xPix, yPix+y, width, alphas, width*y);
			}
			data.alphaData = alphas;
		}
	} else {
		// NOTE Assumes 24 Bit Images
		final int[] pixels = new int[width];
		data = new ImageData(width, height, imageData.depth, new PaletteData(0xff0000, 0x00ff00, 0x0000ff));
		for (int y = 0; y < height; y++) {					
			imageData.getPixels(xPix, y+yPix, width, pixels, 0);
			data.setPixels(0, y, width, pixels, 0);
		}
		if (imageData.alphaData != null) {
			final byte[] alphas = new byte[width];
			for (int y = 0; y < height; y++) {					
				imageData.getAlphas(xPix, y+yPix, width, alphas, 0);
				data.setAlphas(0, y, width, alphas, 0);
			}
		}
	}
	data.alpha = imageData.alpha;
	return data;
}

private final double[] getImageCoords(int bin, boolean pixels, double[] x, double[] y) {
	
	Range xRange = xAxis.getRange();
	Range yRange = yAxis.getRange();

	double minX = xRange.getLower()/bin;
	double minY = yRange.getLower()/bin;
	double maxX = xRange.getUpper()/bin;
	double maxY = yRange.getUpper()/bin;

	// Make sure we have the min and max right
	if(maxX < minX){
		double temp = maxX;
		maxX = minX;
		minX = temp;
	}
	if(maxY < minY){
		double temp = maxY;
		maxY = minY;
		minY = temp;
	}
			
	// Bind the extent of the images to the actual data
	double[] minMaxX= getMinMax(x);
	double[] minMaxY= getMinMax(x);
	
	double minXData = minMaxX[0]/bin;
	minX = Math.max(minXData, minX);
	
	double maxXData = minMaxX[1]/bin;
	maxX = Math.min(maxXData, maxX);
	
	double minYData = minMaxY[0]/bin;
	minY = Math.max(minYData, minY);

	double maxYData = minMaxY[1]/bin;
	maxY = Math.min(maxYData, maxY);

	if (pixels) {
		int x1 = xAxis.getValuePosition(minX, false);
		int y1 = yAxis.getValuePosition(minY, false);
		int x2 = xAxis.getValuePosition(maxX, false);
		int y2 = yAxis.getValuePosition(maxY, false);
		return new double[]{x1, y1, x2, y2};
		
	} else {
	    return new double[]{minX, minY, maxX, maxY};
	}

}

	private double[] getMinMax(double[] data) {
		double min = Double.POSITIVE_INFINITY;
		double max = Double.NEGATIVE_INFINITY;
		
		for (int i = 0; i < data.length; i++) {
			double val = data[i];
			if (val > max) max = val;
			if (val < min) min = val;
		}
		
		return new double[] {min,max};
	}

}
