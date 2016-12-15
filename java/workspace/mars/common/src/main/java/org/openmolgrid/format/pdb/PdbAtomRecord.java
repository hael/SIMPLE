package org.openmolgrid.format.pdb;

import org.apache.log4j.Logger;

/**
 * Class represents a ATOM record format in a pdb file
 * 
 * <table border="4">
 * <tbody>
 * <tr align="top">
 * <th>Name</th>
 * <th>Start</th>
 * <th>End</th>
 * <th>Format</th>
 * <th>Description</th>
 * </tr>
 * <tr align="top">
 * <td>recname</td>
 * <td>0</td>
 * <td>5</td>
 * <td>A6</td>
 * <td>a literal <tt>"ATOM&nbsp;&nbsp;"</tt> (note two trailing spaces).</td>
 * </tr>
 * <tr align="top">
 * <td>serial</td>
 * <td>6</td>
 * <td>10</td>
 * <td>I5</td>
 * <td>atom serial number, e.g. <tt>"&nbsp;&nbsp;&nbsp;86"</tt>. <a href="#pdb-atom-index-notes">See below</a> for
 * details.</td>
 * </tr>
 * <tr align="top">
 * <td></td>
 * <td>11</td>
 * <td>11</td>
 * <td>1X</td>
 * <td>space</td>
 * </tr>
 * <tr align="top">
 * <td>atom</td>
 * <td>12</td>
 * <td>15</td>
 * <td>A4</td>
 * <td>Atom role name, e.g. <tt>"&nbsp;CG1;"</tt>. <a href="#pdb-atom-name-anomalies">See below</a> for details.</td>
 * </tr>
 * <tr align="top">
 * <td>altLoc</td>
 * <td>16</td>
 * <td>16</td>
 * <td>A1</td>
 * <td><a name="atom-variant"> atom variant, officially called the "alternate location indicator".  This is usually
 * <tt>"&nbsp;"</tt> for atoms with well-defined positions, as in this case, but sometimes "A", "B", etc. </a><a
 * href="#pdb-atom-variant-anomalies">See below</a> for details.</td>
 * </tr>
 * <tr align="top">
 * <td>resName</td>
 * <td>17</td>
 * <td>19</td>
 * <td>A3</td>
 * <td>amino acid abbreviation, e.g. <tt>"ARG"</tt>. <a href="#pdb-aa-abbrev-notes">See below</a> for details.</td>
 * </tr>
 * <tr align="top">
 * <td></td>
 * <td>20</td>
 * <td>20</td>
 * <td>1X</td>
 * <td>space</td>
 * </tr>
 * <tr align="top">
 * <td>chainID</td>
 * <td>21</td>
 * <td>21</td>
 * <td>A1</td>
 * <td><a name="pdb-chain-id"> chain ID, usually <tt>"&nbsp;"</tt>, but often <tt>"A"</tt>, <tt>"B"</tt>, etc, for
 * multichain entries. </a><a href="#pdb-atom-chain-notes">See below</a> for details.</td>
 * </tr>
 * <tr align="top">
 * <td>Seqno</td>
 * <td>22</td>
 * <td>26</td>
 * <td>A5</td>
 * <td>residue sequence number (I4) and insertion code (A1), e.g. "&nbsp;&nbsp;11&nbsp;" or "&nbsp;256C". <a
 * href="#pdb-sequence-number">See below</a> for details.</td>
 * </tr>
 * <tr align="top">
 * <td></td>
 * <td>27</td>
 * <td>29</td>
 * <td>3X</td>
 * <td>three spaces</td>
 * </tr>
 * <tr align="top">
 * <td>x</td>
 * <td>30</td>
 * <td>37</td>
 * <td>F8.3</td>
 * <td>atom X coordinate</td>
 * </tr>
 * <tr align="top">
 * <td>y</td>
 * <td>38</td>
 * <td>45</td>
 * <td>F8.3</td>
 * <td>atom Y coordinate</td>
 * </tr>
 * <tr align="top">
 * <td>z</td>
 * <td>46</td>
 * <td>53</td>
 * <td>F8.3</td>
 * <td>atom Z coordinate</td>
 * </tr>
 * <tr align="top">
 * <td>occupancy</td>
 * <td>54</td>
 * <td>59</td>
 * <td>F6.2</td>
 * <td><a name="atom-multiplicity"> atom occupancy, usually <tt>"&nbsp;&nbsp;1.00"</tt>. The sum of </a><a
 * href="#pdb-atom-occupancy-notes">atom occupancies</a> for all <a href="#atom-variant">variants in field&nbsp;4</a>
 * generally add to 1.0.</td>
 * </tr>
 * <tr align="top">
 * <td>tempFactor</td>
 * <td>60</td>
 * <td>65</td>
 * <td>F6.2</td>
 * <td>B value or temperature factor, e.g. <tt>"&nbsp;17.72"</tt>. (I don't use this value, so have nothing to add; see
 * <a href="http://www.rcsb.org/pdb/docs/format/pdbguide2.2/part_62.html">the ATOM record specification</a> discussion
 * of B factors, etc. -- rgr, 8-Oct-96.)</td>
 * </tr>
 * <tr align="top">
 * <td></td>
 * <td>66</td>
 * <td>71</td>
 * <td>A6</td>
 * <td>ignored. [Some older PDB files have footnote numbers here, but this field is not described in the Format 2.1
 * specification. -- rgr, 22-Jan-99.]</td>
 * </tr>
 * <tr align="top">
 * <td>recID</td>
 * <td>72</td>
 * <td>79</td>
 * <td>A8</td>
 * <td><font color="green">[prior to format version 2.0.]</font> record identification field, e.g.
 * <tt>"1AAK&nbsp;146"</tt> (tres FORTRAN, n'est-ce pas?).</td>
 * </tr>
 * <tr align="top">
 * <td>segID</td>
 * <td>72</td>
 * <td>75</td>
 * <td>A4</td>
 * <td>segment identifier, left-justified. <font color="green">[format version 2.0 and later.]</font></td>
 * </tr>
 * <tr align="top">
 * <td>element</td>
 * <td>76</td>
 * <td>77</td>
 * <td>A2</td>
 * <td>element symbol, right-justified. <font color="green">[format version 2.0 and later.]</font></td>
 * </tr>
 * <tr align="top">
 * <td>charge</td>
 * <td>78</td>
 * <td>79</td>
 * <td>A2</td>
 * <td>charge on the atom. <font color="green">[format version 2.0 and later.]</font></td>
 * </tr>
 * </tbody>
 * </table>
 * 
 * 
 * @see http://deposit.rcsb.org/adit/docs/pdb_atom_format.html
 * @author Stefan Bozic
 */
public class PdbAtomRecord {

	/** Logger */
	private static Logger log = Logger.getLogger(PdbAtomRecord.class.getName());

	String recname = null;

	String serial = null;
	String atom = null;
	String altLocation = null;
	String resName = null;
	String chainId = null;
	String seqNo = null;
	String x, y, z = null;
	String occupancy = null;
	String tempFactor = null;
	String recId = null;
	String segId = null;
	String element = null;
	String charge = null;

	/**
	 * @return the recname
	 */
	public String getRecname() {
		return recname;
	}

	/**
	 * @return the serial
	 */
	public String getSerial() {
		return serial;
	}

	/**
	 * @return the atom
	 */
	public String getAtom() {
		return atom;
	}

	/**
	 * @return the altLocation
	 */
	public String getAltLocation() {
		return altLocation;
	}

	/**
	 * @return the resName
	 */
	public String getResName() {
		return resName;
	}

	/**
	 * @return the chainId
	 */
	public String getChainId() {
		return chainId;
	}

	/**
	 * @return the seqNo
	 */
	public String getSeqNo() {
		return seqNo;
	}

	/**
	 * @return the x
	 */
	public String getX() {
		return x;
	}

	/**
	 * @return the y
	 */
	public String getY() {
		return y;
	}

	/**
	 * @return the z
	 */
	public String getZ() {
		return z;
	}

	/**
	 * @return the occupancy
	 */
	public String getOccupancy() {
		return occupancy;
	}

	/**
	 * @return the tempFactor
	 */
	public String getTempFactor() {
		return tempFactor;
	}

	/**
	 * @return the recId
	 */
	public String getRecId() {
		return recId;
	}

	/**
	 * @return the segId
	 */
	public String getSegId() {
		return segId;
	}

	/**
	 * @return the element
	 */
	public String getElement() {
		return element;
	}

	/**
	 * @return the charge
	 */
	public String getCharge() {
		return charge;
	}

	/**
	 * Constructor.
	 * 
	 * @param record
	 */
	public PdbAtomRecord(String record) {
		try {
			recname = record.substring(0, 6).trim();
			serial = record.substring(6, 11).trim();
			atom = record.substring(12, 16).trim();
			altLocation = record.substring(16, 17).trim();
			resName = record.substring(17, 20).trim();
			chainId = record.substring(21, 22).trim();
			seqNo = record.substring(22, 27).trim();
			x = record.substring(30, 38).trim();
			y = record.substring(38, 46).trim();
			z = record.substring(46, 54).trim();
			occupancy = record.substring(54, 60).trim();
			tempFactor = record.substring(60, 66).trim();
			segId = record.substring(72, 76).trim();
			element = record.substring(76, 78).trim();
			charge = record.substring(78, 80).trim();
		} catch (IndexOutOfBoundsException e) {
			log.debug("Record length is too small: " + record.length()
					+ "! It should be 80. Some properties could not be read.", e);
		}
	}
}
