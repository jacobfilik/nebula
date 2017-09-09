package org.eclipse.nebular.visualisation.xygraph.plotlib;

import java.util.Arrays;

import org.eclipse.draw2d.Graphics;
import org.eclipse.draw2d.SWTGraphics;
import org.eclipse.draw2d.geometry.Rectangle;
import org.eclipse.nebula.visualization.xygraph.figures.Axis;
import org.eclipse.nebula.visualization.xygraph.linearscale.Range;
import org.eclipse.swt.SWT;
import org.eclipse.swt.graphics.Color;
import org.eclipse.swt.graphics.Image;
import org.eclipse.swt.graphics.ImageData;
import org.eclipse.swt.graphics.PaletteData;
import org.eclipse.swt.widgets.Display;

public class AxesImage extends AbstractPlotArtist {

	private double[] image;
	private int[] shape;
	double[] xExtent;
	double[] yExtent;
	private double maxMag = 1;
	private double minMag = 0;
	private int currentDownSampleBin = 1;
	private Axis xAxis;
	private Axis yAxis;
	private PaletteData palette = new PaletteData(0xff, 0xff00, 0xff0000);
	private ScaledImageData scaledData;
	
	private Image testImage;

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
		
//			boolean draw = buildImageRelativeToAxes(null);
			if (true) {
				Rectangle bounds2 = axes.getPlotArea().getBounds();			
				if (graphics instanceof SWTGraphics) {
					((SWTGraphics)graphics).setInterpolation(SWT.NONE);
				}
				Rectangle[] sD = calculateSrcDest();
				
				if (sD != null) graphics.drawImage(testImage, sD[0], sD[1]);
//				graphics.drawImage(scaledData.getScaledImage(), scaledData.getXPosition(), scaledData.getYPosition());
			}
		
		
		graphics.popState();
	}
	
	private Rectangle[] calculateSrcDest() {
		Range xRange = xAxis.getLocalRange();
		Range yRange = yAxis.getLocalRange();
		
		double xAxValLow = xRange.getLower();
		double xAxValUp = xRange.getUpper();
		double yAxValLow = yRange.getLower();
		double yAxValUp = yRange.getUpper();
		
		double xAxPixLow = xAxis.getValuePrecisePosition(xAxValLow,false);
		double xAxPixUp= xAxis.getValuePrecisePosition(xAxValUp,false);
		double yAxPixLow = yAxis.getValuePrecisePosition(yAxValLow,false);
		double yAxPixUp = yAxis.getValuePrecisePosition(yAxValUp,false);
		
		
		
		double[] xAxEx = new double[] {xAxValLow,xAxValUp};
		sort(xAxEx);
		
		double[] yAxEx = new double[] {yAxValLow,yAxValUp};
		sort(yAxEx);
		
		double[] xAxPixEx = new double[] {xAxPixLow,xAxPixUp};
		sort(xAxPixEx);
		
		double[] yAxPixEx = new double[] {yAxPixLow,yAxPixUp};
		sort(yAxPixEx);
		
		double[] xEx = xExtent.clone();
		sort(xEx);

		double[] yEx = yExtent.clone();
		sort(yEx);
		
		double dDataPixX = (xEx[1]-xEx[0])/shape[1];
		double dDataPixY = (yEx[1]-yEx[0])/shape[0];
		
		if (xAxEx[0] > xEx[1] || xAxEx[1] < xEx[0] || yAxEx[0] > yEx[1] || yAxEx[1] < yEx[0]) {
			return null;
		}
		
		//val per pixel
		double dx = (xAxEx[1]-xAxEx[0])/(xAxPixEx[1]-xAxPixEx[0]);
		double dy = (yAxEx[1]-yAxEx[0])/(yAxPixEx[1]-yAxPixEx[0]);
		
		int srcX0 = 0;
		int srcY0 = 0;
		int srcX1 = shape[1];
		int srcY1 = shape[0];
		
		int destX0 = (int)(xAxPixLow + (xEx[0]-xAxValLow)/dx);
		int destY0 = (int)(yAxPixUp + (yAxValUp-yEx[1])/dy);
		int destX1 = (int) (xAxPixUp - (xAxValUp-xEx[1])/dx);
		int destY1 = (int) (yAxPixLow - (yEx[0]-yAxValLow)/dy);		
		
		double dx0 = xEx[0]-xAxEx[0];
		double dx1 = xEx[1]-xAxEx[1];
		
		double dy1 = yEx[0]-yAxEx[0];
		double dy0 = yEx[1]-yAxEx[1];
		
		if (dx0 < 0) {
			srcX0 = -1*(int)(dx0/dDataPixX);
			destX0 += (srcX0*dDataPixX/dx);
			if (srcX0 > srcX1) return null;
		}
		
		if (dx1 > 0) {
			int dsrcX1 = -1*(int)(dx1/dDataPixX);
			srcX1 = srcX1 + dsrcX1;
			destX1 += (dsrcX1*dDataPixX/dx);
		}
		
		if (dy0 > 0) {
			srcY0 = (int)(dy0/dDataPixY);
			destY0 += (srcY0*dDataPixY/dy);
		
		}
		
		if (dy1 < 0) {
			int dsrcY1 = (int)(dy1/dDataPixY);
			srcY1 = srcY1 + dsrcY1;
			destY1 += (dsrcY1*dDataPixY/dy);
		}
		
		Rectangle dest = new Rectangle(destX0, 
				destY0, 
				destX1-destX0, 
				destY1-destY0);

		Rectangle src1 = new Rectangle(srcX0, srcY0, srcX1-srcX0, srcY1-srcY0);
//		System.out.println("dest :" + dest.toString());
//		System.out.println("src :" +src1.toString());
		return new Rectangle[] {src1, dest};
		
	}
	
	private void sort(double[] array) {
		if (array[0] > array[1]) {
			double tmp = array[0];
			array[0] = array[1];
			array[1] = tmp;
		}
	}

	/**
	 *
	 * @param xloc  - location of center of arrow
	 * @param yloc  - location of center of arrow
	 * @param length  - size of arrow, normalized from magnitude of vector.
	 * @param theta - anti-clockwise angle from 12 O'clock in radians.
	 */
		
	public void setData(double[] image, int[] shape, double[] extent) {
		this.image = image;
		this.shape = shape;
		if (extent == null) {
			xExtent = new double[]{0, shape[1]};
			yExtent = new double[]{0, shape[0]};
		} else {
			xExtent = new double[]{extent[0], extent[1]};
			yExtent = new double[]{extent[2], extent[3]};
		}
		
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
		
//		ImageData imageData = new ImageData(shape[1], shape[0], 8, palette, 1, buffer);
//		new ImageData(width, height, depth, palette, scanlinePad, data)
		
		makeTestImage();
		buildImageRelativeToAxes(null);
	}



	@Override
	protected void onAxesSet() {
		this.xAxis = axes.getPrimaryXAxis();
		this.yAxis = axes.getPrimaryYAxis();
		buildImageRelativeToAxes(null);
	}
	
	private void makeTestImage() {
		
		byte[] buffer = new byte[image.length];
		
		for (int i = 0; i < image.length; i++) {
			buffer[i] = (byte)(((image[i]-minMag)/(maxMag-minMag))*255);
		}
		
		ImageData data = new ImageData(shape[1], shape[0], 8, palette, 1, buffer);
		
		testImage = new Image(Display.getDefault(), data);
		
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
		
		ImageData data = new ImageData(shape[1], shape[0], 8, palette, 1, buffer);
		
		if (xAxis == null || yAxis == null) return false;
		
		
		Range xRange = xAxis.getRange();
		Range yRange = yAxis.getRange();
		
		boolean noSlice = noSliceNeeded(xRange, yRange, xExtent, yExtent);
		
		if (noSlice) {
			System.out.println("no slice needed");
		}
		
		int x0 = xAxis.getValuePosition(xExtent[0], false);
		int x1 = xAxis.getValuePosition(xExtent[1], false);
		int y0 = yAxis.getValuePosition(yExtent[0], false);
		int y1 = yAxis.getValuePosition(yExtent[1], false);
		
		boolean xAxisInc = xRange.getLower() < xRange.getUpper();
		boolean yAxisInc = yRange.getLower() < yRange.getUpper();
		
		boolean xDataInc = xExtent[0] < xExtent[1];
		boolean yDataInc = xExtent[0] < xExtent[1];
		
		double xImagePixelSize = (xExtent[1]-xExtent[0])/shape[1];
		double yImagePixelSize = (yExtent[1]-yExtent[0])/shape[0];
		
		//calculate slicing, should be 1 pixel larger than plot area for smooth scrolling
		//Slice x and y starts
		int flooredYStart = (int)Math.floor((yRange.getLower()-yExtent[0])/yImagePixelSize);
		int flooredXStart = (int)Math.floor((xRange.getLower()-xExtent[0])/xImagePixelSize);
		System.out.println(flooredXStart);
		
		//slice stops - todo
		
		
		//Get the axes coodinates visible on screen
		double[] da = getImageCoords(1, false, xExtent, yExtent);
		double minX = da[0];
		double minY = da[1];
		double maxX = da[2];
		double maxY = da[3];
		
		//Get the data points per axis step
//		double xAxValPerPoint =  Math.abs(getAxisValuePerDataPoint(minX,maxX,xData));
//		double yAxValPerPoint =  Math.abs(getAxisValuePerDataPoint(minY,maxY,yData));
//		
//		if (Double.isNaN(xAxValPerPoint)|| Double.isNaN(yAxValPerPoint)) {
//			return false;
//		}
//		
//		//get x and y start position in data array (floored)
//		int xPix = getPositionInAxis(xDataInc ? minX : maxX,xData,true)/currentDownSampleBin;
//		int yPix = getPositionInAxis(yDataInc ? minY : maxY,yData,true)/currentDownSampleBin;
//		
//		//determine sub pixel offset in data frame
//		double xdiff = (xData[xPix])-(xDataInc ? minX : maxX);
//		double ydiff = (yData[yPix])-(yDataInc ? minY : maxY);
//		
//		//Determine the corresponding number of data points in x and y (floor min, ceil max)
//		int xDataPoints = getNumberOfDataPoints(minX,maxX,xData, xDataInc);
//		int yDataPoints = getNumberOfDataPoints(minY,maxY,yData, yDataInc);
//		
//		int xDataPointsDS = xDataPoints/currentDownSampleBin;
//		int yDataPointsDS = yDataPoints/currentDownSampleBin;
//		
//		//Get matching screen co-ordinates in pixels
//		double[] screenCoords = getImageCoords(1, true, xData, yData);
		
		//calculate full number of data points
//		double realx = (maxX-minX)/xAxValPerPoint;
//		double realy = (maxY-minY)/yAxValPerPoint;
//		
//		// Ratio of screen pixels to full number of data points
//		double xScale = Math.abs((screenCoords[2]-screenCoords[0]) / (realx));
//		double yScale = Math.abs((screenCoords[3]-screenCoords[1]) / (realy));
		
		//FIXME Handle minimum cases
		
		// Size of image data in index
//		final int xSize = imageData.width;
//		final int ySize = imageData.height;
		
		double[] minMaxX = xExtent.clone();
		Arrays.sort(minMaxX);
		double[] minMaxY = xExtent.clone();
		Arrays.sort(minMaxY);
		
		double xmin = minMaxX[0];
		double xmax = minMaxX[1];
		double ymin = minMaxY[0];
		double ymax = minMaxY[1];
		
		double xLen = shape[1];
		double yLen = shape[0];
		
		double xp = ((xmax-xmin)/(xLen-1))/2;
		double yp = ((ymax-ymin)/(yLen-1))/2;
		
//		double xOffset = (((-xdiff+xp)/xAxValPerPoint))*xScale;
//		double yOffset = (((-ydiff+yp)/yAxValPerPoint))*yScale;
////		xOffset = 0;
//		
//
//		//FIXME origins
//		
//		int width  = xPix+xDataPointsDS > xSize  ? Math.min(xSize, xDataPointsDS) : xDataPointsDS;
//		int height = yPix+yDataPointsDS > ySize  ? Math.min(ySize, yDataPointsDS) : yDataPointsDS;
//		
//		int scaleWidth  = Math.max(1, (int) (xDataPoints*xScale));
//		int scaleHeight = Math.max(1, (int) (yDataPoints*yScale));
		
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
//			data = sliceImageData(data, width, height, xPix, yPix, ySize);

			// create the scaled image
			// We are suspicious if the algorithm wants to create an image
			// bigger than the screen size and in that case do not scale
			// Fix to http://jira.diamond.ac.uk/browse/SCI-926
			boolean proceedWithScale = true;
			try {
//				if (screenRectangle == null) {
//					screenRectangle = Display.getCurrent().getPrimaryMonitor().getClientArea();
//				}
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
			
//			if (!xDataInc) scaleWidth*=-1;
//			if (!yDataInc) scaleHeight*=-1;
			
			int scaleWidth = (x1-x0);
			int scaleHeight = (y1-y0);
			
			if (scaleWidth == 0 || scaleHeight == 0) return false;
			
			Image scaledImage = null;
			if (proceedWithScale) {
				data = data!=null ? data.scaledTo(scaleWidth, scaleHeight) : null;
				scaledImage = data!=null ? new Image(Display.getDefault(), data) : null;
			} else if (scaledImage==null) {
				scaledImage = data!=null ? new Image(Display.getDefault(), data) : null;
			}
			
			scaledData.setX(x0);
			scaledData.setY(y1);
			scaledData.setXoffset(0);
			scaledData.setYoffset(0);
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
	
	private boolean noSliceNeeded(Range x, Range y, double[] xExtent, double[] yExtent) {
		
		double[] xAxis = new double[] {x.getLower(), x.getUpper()};
		Arrays.sort(xAxis);
		
		double[] yAxis = new double[] {y.getLower(), y.getUpper()};
		Arrays.sort(yAxis);
		
		double[] yExt = yExtent.clone();
		Arrays.sort(yExt);
		
		double[] xExt = xExtent.clone();
		Arrays.sort(xExt);
		
		return yAxis[0] <= yExt[0] && yAxis[1] >= yExt[1] && xAxis[0] <= xExt[0] && xAxis[1] >= xExt[1]; 
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

private final double[] getImageCoords(int bin, boolean pixels, double[] xExtent, double[] yExtent) {
	
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
	double[] minMaxX= xExtent[0] < xExtent[1] ? xExtent.clone() : new double[]{xExtent[1], xExtent[0]};
	double[] minMaxY= yExtent[0] < yExtent[1] ? yExtent.clone() : new double[]{yExtent[1], yExtent[0]};
	
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
