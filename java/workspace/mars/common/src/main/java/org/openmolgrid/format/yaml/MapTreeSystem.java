/**
 * 
 */
package org.openmolgrid.format.yaml;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import javax.swing.event.TreeModelListener;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreeModel;
import javax.swing.tree.TreePath;

/**
 * @author Frederic Bonnet
 *
 */
public class MapTreeSystem { //implements TreeModel {

	private DefaultMutableTreeNode root;
	private DefaultMutableTreeNode node;
	private DefaultMutableTreeNode node1;
	private DefaultMutableTreeNode node2;
	private DefaultMutableTreeNode node3;
	//private DefaultMutableTreeNode parentNode;
	private DefaultMutableTreeNode childrenNode;
	private Map<?, ?> mapNextLevel;
	private Vector listeners = new Vector();
	private Set<?> set = null;
	private Object[] arrayNextLevel;
	//starts of the methods
	public MapTreeSystem(DefaultMutableTreeNode rootsys, Map<?, ?> mapNextLevel) {
		root = rootsys;
		this.mapNextLevel = mapNextLevel;
		set = mapNextLevel.keySet();
	}
	private Object[] getObject(Set<?> setin) {
		arrayNextLevel = new String[setin.size()];
		arrayNextLevel = setin.toArray();
		return arrayNextLevel;
	}
	public boolean isValideMap(Object[] object, int elemNextLevel) {
		boolean result = false;
		if (mapNextLevel.get(object[elemNextLevel]) instanceof Map<?, ?> ) {
//			System.out.println("Map<?, ?>"+mapNextLevel.get(object[elemNextLevel]));
			result = true;
		} else if (mapNextLevel.get(object[elemNextLevel]) instanceof List<?> ) {
//			System.out.println("List<?>"+mapNextLevel.get(object[elemNextLevel]));
			result = false;
		} else if (mapNextLevel.get(object[elemNextLevel]) instanceof String ) {
//			System.out.println("String"+mapNextLevel.get(object[elemNextLevel]));
			result = false;
		}
		return result;
	}
    public Object getRoot() {
    	return root;
    }
    public void setRoot(DefaultMutableTreeNode root) {
    	this.root = root;
    }
    public Object[] getFirstChild(Object[] parentObject, int index, int elemNextLevel0) {
		if ( !isValideMap(parentObject, elemNextLevel0) ) {
			throw new IllegalArgumentException("parent is not a node in this model");
		}
		if (mapNextLevel.get(parentObject[elemNextLevel0]) instanceof Map<?, ?> ) {
			Map<?, ?> firstChildrenMap = (Map<?, ?>)mapNextLevel.get(parentObject[elemNextLevel0]);
			if ( firstChildrenMap != null ) {
				childrenNode /**node*/ = new DefaultMutableTreeNode(parentObject[elemNextLevel0].toString());
				root.add(childrenNode);
				//this.parentNode = root;
				
				Set<?> firstChildrenSet = firstChildrenMap.keySet();
        		for ( int elemNextLevel1 = 0 ; elemNextLevel1 < firstChildrenSet.size() ; elemNextLevel1++ ) {
        			getChild(childrenNode, firstChildrenMap, getObject(firstChildrenSet), index, elemNextLevel1);
    				index++;
        		}
			}
		} else if (mapNextLevel.get(parentObject[elemNextLevel0]) instanceof List<?> ) {
			List<?> firstChildrenList = (List<?>)mapNextLevel.get(parentObject[elemNextLevel0]);
			if ( firstChildrenList != null ) {
				childrenNode /**node*/ = new DefaultMutableTreeNode(firstChildrenList);
				root.add(childrenNode);
				//this.parentNode = root;
			}
		} else if (mapNextLevel.get(parentObject[elemNextLevel0]) instanceof String ) {
			String firstChildrenString = (String)mapNextLevel.get(parentObject[elemNextLevel0]);
			if ( firstChildrenString != null ) {
				childrenNode /**node*/ = new DefaultMutableTreeNode(firstChildrenString);
				root.add(childrenNode);
				//this.parentNode = root;
			}
		} else if (mapNextLevel.get(parentObject[elemNextLevel0]) instanceof Integer ) {
			Integer firstChildrenInteger = (Integer)mapNextLevel.get(parentObject[elemNextLevel0]);
			if ( firstChildrenInteger != null ) {
				childrenNode /**node*/ = new DefaultMutableTreeNode(firstChildrenInteger);
				root.add(childrenNode);
				//this.parentNode = root;
			}
		}		
		return arrayNextLevel;
	}
	public Object getChild(DefaultMutableTreeNode parentNodeIn , Map<?, ?> parentMap, Object[] parent, int index, int elemNextLevel) {
		node = new DefaultMutableTreeNode(parent[elemNextLevel].toString());
		parentNodeIn.add(node);
		if (parentMap.get(parent[elemNextLevel]) instanceof Map<?, ?> ) {
			Map<?, ?> childrenMap = (Map<?, ?>)parentMap.get(parent[elemNextLevel]);
			if ( childrenMap != null ) {
				Set<?> set = childrenMap.keySet();
        		for ( int elemNextLevel1 = 0 ; elemNextLevel1 < set.size() ; elemNextLevel1++ ) {
        			node1 = new DefaultMutableTreeNode(getObject(set)[elemNextLevel1].toString());
        			node.add(node1);
//        			getChild(node1, childrenMap, getObject(set), index, elemNextLevel1);
        			if (childrenMap.get(getObject(set)[elemNextLevel1]) instanceof Map<?, ?>) {
        				Map<?, ?> grandChildrenMap = (Map<?, ?>)childrenMap.get(getObject(set)[elemNextLevel1]);
        				Set<?> set1 = grandChildrenMap.keySet();
        				if (grandChildrenMap != null ) {
        	        		for ( int elemNextLevel2 = 0 ; elemNextLevel2 < set1.size() ; elemNextLevel2++ ) {
        	        			node2 = new DefaultMutableTreeNode(getObject(set1)[elemNextLevel2].toString());
        	        			node1.add(node2);
        	        			getChild(grandChildrenMap,set1,elemNextLevel2);
        	        		}
        				}
        			} else if (childrenMap.get(getObject(set)[elemNextLevel1]) instanceof List<?>) {
        				List<?> grandChildrenList = (List<?>)childrenMap.get(getObject(set)[elemNextLevel1]);
        				if (grandChildrenList != null ) {
    	        			node2 = new DefaultMutableTreeNode(grandChildrenList);
    	        			node1.add(node2);
        				}
        			} else if (childrenMap.get(getObject(set)[elemNextLevel1]) instanceof String) {
        				String grandChildrenString = (String)childrenMap.get(getObject(set)[elemNextLevel1]);
        				if (grandChildrenString != null ) {
    	        			node2 = new DefaultMutableTreeNode(grandChildrenString.toLowerCase());
    	        			node1.add(node2);
        				}
        			} else if (childrenMap.get(getObject(set)[elemNextLevel1]) instanceof Integer) {
        				Integer grandChildrenInteger = (Integer)childrenMap.get(getObject(set)[elemNextLevel1]);
        				if (grandChildrenInteger != null ) {
    	        			node2 = new DefaultMutableTreeNode(grandChildrenInteger);
    	        			node1.add(node2);
        				}
        			}
        			//
        		}
			}
		} else if (parentMap.get(parent[elemNextLevel]) instanceof List<?>) {
			List<?> childrenList = (List<?>)parentMap.get(parent[elemNextLevel]);
			if ( childrenList != null ) {
				node1/**node*/ = new DefaultMutableTreeNode(childrenList.toString());
				node.add(node1);
			}
		} else if (parentMap.get(parent[elemNextLevel]) instanceof String) {
			String childrenString = (String)parentMap.get(parent[elemNextLevel]);
			if ( childrenString != null ) {
				node1/**node*/ = new DefaultMutableTreeNode(childrenString.toLowerCase());
				node.add(node1);
			}
		} else if (parentMap.get(parent[elemNextLevel]) instanceof Integer ) {
			Integer childrenInteger = (Integer)parentMap.get(parent[elemNextLevel]);
			if ( childrenInteger != null ) {
				node1/**node*/ = new DefaultMutableTreeNode(childrenInteger);
				node.add(node1);
			}
		}
		return null;
	}
	/**
	 * here the recursion is not the perfect with the node allocation 
	 */
	public Object getChild(Map<?, ?> childrenMap, Set<?> set1, int elemNextLevel2) {
		//TODO: must fix the node allocation in the recursion mismatch in the node.
		if (childrenMap.get(getObject(set1)[elemNextLevel2]) instanceof Map<?, ?>) {
			Map<?, ?> grandChildrenMap = (Map<?, ?>)childrenMap.get(getObject(set1)[elemNextLevel2]);
			Set<?> set2 = grandChildrenMap.keySet();
			if (grandChildrenMap != null ) {
        		for ( int elemNextLevel3 = 0 ; elemNextLevel3 < set2.size() ; elemNextLevel3++ ) {
        			node3 = new DefaultMutableTreeNode(getObject(set2)[elemNextLevel3].toString());
        			node2.add(node3);
        			getChild(grandChildrenMap, set2, elemNextLevel3);
        		}
			}
		} else if (childrenMap.get(getObject(set1)[elemNextLevel2]) instanceof List<?>) {
			List<?> grandChildrenList = (List<?>)childrenMap.get(getObject(set1)[elemNextLevel2]);
			if (grandChildrenList != null ) {
    			node3 = new DefaultMutableTreeNode(grandChildrenList);
    			node2.add(node3);
			}
		} else if (childrenMap.get(getObject(set1)[elemNextLevel2]) instanceof String) {
			String grandChildrenString = (String)childrenMap.get(getObject(set1)[elemNextLevel2]);
			if (grandChildrenString != null ) {
    			node3 = new DefaultMutableTreeNode(grandChildrenString.toLowerCase());
    			node2.add(node3);
			}
		} else if (childrenMap.get(getObject(set1)[elemNextLevel2]) instanceof Integer) {
			Integer grandChildrenInteger = (Integer)childrenMap.get(getObject(set1)[elemNextLevel2]);
			if (grandChildrenInteger != null ) {
    			node3 = new DefaultMutableTreeNode(grandChildrenInteger);
    			node2.add(node3);
			}
		}
		return null;
	}
	public int getIndexOfChild(Object parent , Object child) {
		return 0;
	}
    public boolean isLeaf(Object node) {
    	return false;
    }
	public int getChildCount(Object parent) {
		return 0;
	}
	public void valueForPathChanged(TreePath path, Object value) {
		
	}
	public void addTreeModelListener(TreeModelListener listener) {
		listeners.add(listener);
	}
	public void removeTreeModelListener(TreeModelListener listener) {
		listeners.remove(listener);
	}
}
