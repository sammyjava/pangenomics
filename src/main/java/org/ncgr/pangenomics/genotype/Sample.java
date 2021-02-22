package org.ncgr.pangenomics.genotype;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import java.util.TreeSet;

/**
 * Encapsulates a sample, which is just a name and a label.
 */
public class Sample implements Comparable {

    public String name = null;
    public String label = null;

    /**
     * Construct from a name and label
     */
    public Sample(String name, String label) {
	this.name = name;
	this.label = label;
    }

    /**
     * Construct from a tab-delimited labels.txt or paths.txt line
     *
     * 0      1     2     3     4   5
     * sample label other stuff can follow
     */
    public Sample(String line) {
	String[] fields = line.split("\t");
	this.name = fields[0];
	this.label = fields[1];
    }

    /**
     * Return case status from the label value.
     */
    public boolean isCase() {
	if (label==null) {
	    return false;
	} else {
	    return label.toLowerCase().equals("case");
	}
    }

    /**
     * Return control status from the label value.
     */
    public boolean isControl() {
	if (label==null) {
	    return false;
	} else {
	    return label.toLowerCase().equals("ctrl") || label.toLowerCase().equals("control");
	}
    }

    /**
     * Return true if label==false
     */
    public boolean isUnlabeled() {
	return label==null;
    }

    /**
     * Compare alphabetically on name.
     */
    @Override
    public int compareTo(Object o) {
	Sample that = (Sample) o;
	return this.name.compareTo(that.name);
    }

    /**
     * Samples are equal if they have the same name.
     */
    @Override
    public boolean equals(Object o) {
	Sample that = (Sample) o;
	return this.name.equals(that.name);
    }

    /**
     * Read an ordered set of Samples from a labels.txt file.
     */
    public static TreeSet<Sample> readSamples(File labelsFile) throws FileNotFoundException, IOException {
	TreeSet<Sample> samples = new TreeSet<>();
        String line = null;
        BufferedReader reader = new BufferedReader(new FileReader(labelsFile));
        while ((line=reader.readLine())!=null) {
            if (line.startsWith("#")) continue;      // comment
	    if (line.startsWith("sample")) continue; // header
	    samples.add(new Sample(line));
        }
        reader.close();
	return samples;
    }

    /**
     * @override
     */
    public String toString() {
	return name+"\t"+label;
    }
}
