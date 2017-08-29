package org.eclipse.nebular.visualisation.xygraph.plotlib;

import java.util.ArrayList;
import java.util.List;

import org.eclipse.draw2d.geometry.Rectangle;
import org.eclipse.nebula.visualization.xygraph.figures.XYGraph;
import org.eclipse.nebula.visualization.xygraph.util.XYGraphMediaFactory;

public class Axes extends XYGraph {

	List<AbstractPlotArtist> artists = new ArrayList<>();
	
	public void addArtist(AbstractPlotArtist artist) {
		artist.setAxes(this);
		if (artist instanceof Line2D) {
			Line2D line = (Line2D)artist;
			if (line.getTraceColor() == null) { // Cycle through default colors
				line.setTraceColor(XYGraphMediaFactory.getInstance().getColor(
						DEFAULT_TRACES_COLOR[plotArea.getTraceList().size() % DEFAULT_TRACES_COLOR.length]));
			}
		}
		if (artist instanceof VectorField) {
			VectorField line = (VectorField)artist;
			
				line.setColor(XYGraphMediaFactory.getInstance().getColor(
						DEFAULT_TRACES_COLOR[plotArea.getTraceList().size() % DEFAULT_TRACES_COLOR.length]));
			
		}
		artists.add(artist);
		this.getPlotArea().add(artist);
		revalidate();
		repaint();
	}
	
	@Override
	protected void layout() {
		super.layout();
		final Rectangle clientArea = getClientArea();
		for (AbstractPlotArtist trace : artists) {
			if (trace != null && trace.isVisible())
				// Shrink will make the trace has no intersection with axes,
				// which will make it only repaints the trace area.
				trace.setBounds(clientArea);// .getCopy().shrink(1, 1));
		}

		super.layout();
	}
	
}
