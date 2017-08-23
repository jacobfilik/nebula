package org.eclipse.nebular.visualisation.xygraph.plotlib;

import org.eclipse.draw2d.Figure;
import org.eclipse.nebula.visualization.xygraph.figures.Axis;

public abstract class AbstractPlotArtist extends Figure {

	protected Axes axes;
	
	
	public void setAxes(Axes axes) {
		this.axes = axes;
		onAxesSet();
	}
	
	protected abstract void onAxesSet();
	
	public Axis getXAxis() {
		return axes.getPrimaryXAxis();
	}
	
	public Axis getYAxis() {
		return axes.getPrimaryYAxis();
	}
}
