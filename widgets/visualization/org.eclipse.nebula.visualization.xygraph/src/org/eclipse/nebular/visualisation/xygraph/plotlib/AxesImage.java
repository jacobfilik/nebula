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
import org.eclipse.swt.graphics.Rectangle;
import org.eclipse.swt.widgets.Display;

public class AxesImage extends AbstractPlotArtist {

	private double[] image;
	private double[] xData;
	private double[] yData;
	private double maxMag = 1;
	private double minMag = 0;
	private int currentDownSampleBin = 1;
	private Axis xAxis;
	private Axis yAxis;
	private PaletteData palette = new PaletteData(0xff, 0xff00, 0xff0000);
	private Color color;
	private Rectangle screenRectangle;
	private ScaledImageData scaledData;

	public AxesImage() {
		
		scaledData = new ScaledImageData();

	}

	
	@Override
	protected void paintFigure(Graphics graphics) {
		
		super.paintFigure(graphics);

		/**
		 * This is not actually needed except that when there
		 * are a number of opens of an image, e.g. when moving
		 * around an h5 gallery with arrow keys, it looks smooth 
		 * with this in.
		 */
		if (scaledData.getScaledImage()==null) {
//			boolean imageReady = buildImageRelativeToAxes(im);
		}
		
	
		graphics.pushState();
		
		// Offsets and scaled image are calculated in the createScaledImage method.
		
			boolean draw = buildImageRelativeToAxes(null);
			if (draw) graphics.drawImage(scaledData.getScaledImage(), scaledData.getXPosition(), scaledData.getYPosition());
		
		
		graphics.popState();
	}

	/**
	 *
	 * @param xloc  - location of center of arrow
	 * @param yloc  - location of center of arrow
	 * @param length  - size of arrow, normalized from magnitude of vector.
	 * @param theta - anti-clockwise angle from 12 O'clock in radians.
	 */
		
	public void setData(double[] image, double[] x, double[] y) {
		this.image = image;
		this.xData = x;
		this.yData = y;
		
		maxMag = Double.NEGATIVE_INFINITY;
		minMag = Double.POSITIVE_INFINITY;
		
		for (int i = 0; i < image.length; i++) {
			double val = image[i];
			if (val > maxMag) maxMag = val;
			
			if (val < minMag) minMag = val;
		}
		
		byte[] buffer = new byte[image.length];
		
		for (int i = 0; i < image.length; i++) {
			buffer[i] = (byte)(((image[i]-minMag)/(maxMag-minMag))*255-128);
		}
		
		ImageData imageData = new ImageData(x.length, y.length, 8, palette, x.length, buffer);
//		new ImageData(width, height, depth, palette, scanlinePad, data)
		
		
		buildImageRelativeToAxes(imageData);
	}



	@Override
	protected void onAxesSet() {
		this.xAxis = axes.getPrimaryXAxis();
		this.yAxis = axes.getPrimaryYAxis();
		buildImageRelativeToAxes(null);
	}
	
	private boolean buildImageRelativeToAxes(ImageData imageData) {
		
		// slice data to get current zoom area
		/**     
		 *      x1,y1--------------x2,y2
		 *        |                  |
		 *        |                  |
		 *        |                  |
		 *      x3,y3--------------x4,y4
		 *      
		 */
		
	byte[] buffer = new byte[image.length];
		
		for (int i = 0; i < image.length; i++) {
			buffer[i] = (byte)(((image[i]-minMag)/(maxMag-minMag))*255);
		}
		
		imageData = new ImageData(xData.length, yData.length, 8, palette, 1, buffer);
		
		ImageData data = imageData;
		
		
		if (xAxis == null || yAxis == null) return false;
		
		
		boolean xDataInc = xData[0] < xData[xData.length-1];
		boolean yDataInc = yData[0] < yData[yData.length-1];
		
		//Get the axes coodinates visible on screen
		double[] da = getImageCoords(1, false, xData, yData);
		double minX = da[0];
		double minY = da[1];
		double maxX = da[2];
		double maxY = da[3];
		
		//Get the data points per axis step
		double xAxValPerPoint =  Math.abs(getAxisValuePerDataPoint(minX,maxX,xData));
		double yAxValPerPoint =  Math.abs(getAxisValuePerDataPoint(minY,maxY,yData));
		
		if (Double.isNaN(xAxValPerPoint)|| Double.isNaN(yAxValPerPoint)) {
			return false;
		}
		
		//get x and y start position in data array (floored)
		int xPix = getPositionInAxis(xDataInc ? minX : maxX,xData,true)/currentDownSampleBin;
		int yPix = getPositionInAxis(yDataInc ? minY : maxY,yData,true)/currentDownSampleBin;
		
		//determine sub pixel offset in data frame
		double xdiff = (xData[xPix])-(xDataInc ? minX : maxX);
		double ydiff = (yData[yPix])-(yDataInc ? minY : maxY);
		
		//Determine the corresponding number of data points in x and y (floor min, ceil max)
		int xDataPoints = getNumberOfDataPoints(minX,maxX,xData, xDataInc);
		int yDataPoints = getNumberOfDataPoints(minY,maxY,yData, yDataInc);
		
		int xDataPointsDS = xDataPoints/currentDownSampleBin;
		int yDataPointsDS = yDataPoints/currentDownSampleBin;
		
		//Get matching screen co-ordinates in pixels
		double[] screenCoords = getImageCoords(1, true, xData, yData);
		
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
		
		double[] minMaxX = getMinMax(xData);
		double[] minMaxY = getMinMax(yData);
		
		double xmin = minMaxX[0];
		double xmax = minMaxX[1];
		double ymin = minMaxY[0];
		double ymax = minMaxY[1];
		
		double xLen = xData.length;
		double yLen = yData.length;
		
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
//				if (scaleWidth>screenRectangle.width*2 ||  scaleHeight>screenRectangle.height*2) {
//					logger.error("Image scaling algorithm has malfunctioned and asked for an image bigger than the screen!");
//					logger.debug("scaleWidth="+scaleWidth);
//					logger.debug("scaleHeight="+scaleHeight);
//					proceedWithScale = false;
//				}
			} catch (Throwable ne) {
				proceedWithScale = true;
			}

			
			Range xlocalRange = yAxis.getLocalRange();
			Range ylocalRange = yAxis.getLocalRange();
			
			if (xlocalRange.getLower() > xlocalRange.getUpper()) {
				xDataInc = !xDataInc;
			}
			
			if (ylocalRange.getLower() < ylocalRange.getUpper()) {
				yDataInc = !yDataInc;
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
			
//			logger.error("Image scaling has malfunctioned");
//			logger.debug("Trace name = "+getName());
//			logger.debug("scaleWidth = "+scaleWidth);
//			logger.debug("scaleHeight = "+scaleHeight);
//			logger.debug("width = "+width);
//			logger.debug("height = "+height);
//			logger.debug("xPix = "+xPix);
//			logger.debug("yPix = "+yPix);
			
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
private final int getNumberOfDataPoints(double min, double max, double[] axis, boolean increasing) {
	
	int minp = getPositionInAxis(min,axis,true);
	int maxp = getPositionInAxis(max,axis,false);
	
	return increasing ? maxp-minp+1 : minp-maxp+1;
}
private final int getPositionInAxis(double val, double[] axis, boolean floor) {
	if (axis.length == 1) return 0;
	int pos = 0;
	double dif = 0;
	double test = Double.MAX_VALUE;
	
	for (int i = 0; i < axis.length; i++) {
		double d = axis[i];
		double ad = Math.abs(d-val);
		if (ad < test) {
			test = ad;
			pos = i;
			dif = d-val;
		}
	}
	
	if (floor) {
		if (dif > 0) pos--;
	} else {
		if (dif < 0) pos++;
	}
	
	return pos;
	
}


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
	double[] minMaxY= getMinMax(y);
	
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
