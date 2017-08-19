package org.eclipse.nebula.visualization.xygraph.figures;

import org.eclipse.draw2d.Figure;

public class AbstractPlotDrawable extends Figure {
	
	private IXYGraph xyGraph;
	private Axis xAxis;
	private Axis yAxis;
	private String name;

	public void setXYGraph(IXYGraph graph) {
		this.xyGraph = graph;
	}
	
	/**
	 * @return the xAxis
	 */
	public Axis getXAxis() {
		return xAxis;
	}
	
	/**
	 * @return the xAxis
	 */
	public Axis getYAxis() {
		return yAxis;
	}
	
	/**
	 * @return the name of the trace
	 */
	public String getName() {
		return name;
	}
}
