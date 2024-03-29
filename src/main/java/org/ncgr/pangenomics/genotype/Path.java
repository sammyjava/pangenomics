package org.ncgr.pangenomics.genotype;

import java.util.List;
import java.util.LinkedList;

import org.jgrapht.Graph;
import org.jgrapht.graph.GraphWalk;

/**
 * An extension of GraphWalk to provide an individual path through the graph and methods appropriate for frequented regions.
 * It must implement Comparable, otherwise identical paths from different individuals will be regarded equal.
 *
 * @author Sam Hokin
 */
public class Path extends GraphWalk<Node,Edge> implements Comparable {

    // the Sample associated with this Path
    private Sample sample;

    /**
     * Create a path defined by a List of Nodes and a Sample. weight=1.0.
     */
    public Path(Graph<Node,Edge> graph, List<Node> nodes, Sample sample) {
        super(graph, nodes, 1.0);
        this.sample = sample;
    }

    /**
     * Return true if the two paths have the same sample.
     */
    @Override
    public boolean equals(Object o) {
        Path that = (Path) o;
        return this.sample.equals(that.sample);
    }

    /**
     * Must override hashCode() for Map keys.
     */
    @Override
    public int hashCode() {
	return this.sample.name.hashCode();
    }

    /**
     * Compare based on sample.
     */
    public int compareTo(Object o) {
        Path that = (Path) o;
        return this.sample.compareTo(that.sample);
    }

    /**
     * Return the sample.
     */
    public Sample getSample() {
	return sample;
    }
     
    /**
     * Return the sample.name
     */
    public String getName() {
	return sample.name;
    }
    
    /**
     * Return the sample.label
     */
    public String getLabel() {
	return sample.label;
    }

    /**
     * Return true if this path belongs to a "case" sample.
     */
    public boolean isCase() {
	return sample.isCase();
    }

    /**
     * Return true if this path belongs to a "control" sample.
     */
    public boolean isControl() {
	return sample.isControl();
    }

    /**
     * Return the nodes that this path traverses, in order of traversal.
     */
    public List<Node> getNodes() {
        return getVertexList();
    }

    /**
     * Return the edges that this path follows along, in order of traversal.
     */
    public List<Edge> getEdges() {
        return getEdgeList();
    }

    /**
     * Return a bespoke string representation, using NodeSet.toString() to list the nodes.
     */
    @Override
    public String toString() {
        NodeSet ns = new NodeSet(getNodes());
    	return sample.toString()+"\t"+ns.toString();
    }

    /**
     * Return the subpath inclusively between the two given nodes (empty if one of the nodes is not present in this path).
     * @param nl the "left" node
     * @param nr the "right" node
     * @return the Path inclusively between nl and nr
     */
    public Path subpath(Node nl, Node nr) {
        List<Node> subnodes = new LinkedList<>();
        List<Node> nodeList = getNodes();
        if (nodeList.contains(nl) && nodeList.contains(nr)) {
            if (nl.equals(nr)) {
                subnodes.add(nl);
            } else {
                boolean started = false;
                boolean finished = false;
                for (Node node : nodeList) {
                    if (!started && node.equals(nl)) {
                        started = true;
                        subnodes.add(node);
                    } else if (node.equals(nr) && !finished) {
                        subnodes.add(node);
                        finished = true;
                    } else if (started && !finished) {
                        subnodes.add(node);
                    }
                }
            }
        }
        return new Path(this.graph, subnodes, sample);
    }

    /**
     * Return true if the given Path represents a subpath of this Path.
     * @param path the path to be compared with this one
     * @return true if path is a subpath of this
     */
    public boolean contains(Path path) {
        List<Node> thisNodes = this.getNodes();
        List<Node> pathNodes = path.getNodes();
        boolean match = false;
        for (Node m : pathNodes) {
            if (!match && thisNodes.contains(m)) {
                match = true;
            } else if (match && !thisNodes.contains(m)) {
                match = false;
                break;
            }
        }
        return match;
    }

    /**
     * Return true if this path traverses the given node (identified by id).
     * @param node the node
     * @return true if this path traverses the node
     */
    public boolean traverses(Node node) {
	return getNodes().contains(node);
    }
}
