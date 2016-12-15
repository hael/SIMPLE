/**
 * 
 */
package org.openmolgrid.format.yaml;

import java.awt.BorderLayout;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Enumeration;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.JTree;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreeNode;
import javax.swing.tree.TreePath;

/**
 * @author Frederic Bonnet
 * @date 15th of June 2012.
 *
 */
public class YamlWindowPoper extends JDialog {

	/**
	 * generating the UID for the extended class.
	 */
	private static final long serialVersionUID = 6875704250387622703L;
	/** 
	 * the global variables for the class, needed for the constructor
	 */
	private JTree timeYamlFullTree; 
	private JTree timeYamlRecursivelyTree;
	private DefaultMutableTreeNode treeNodeTimeYamlParse_RecursivelyTree_root;
    private DefaultTreeModel m_modelRecursivelyTree;  

	/**
	 * The private variables used in the method
	 */
	/** variables for the Panel and GUI */
	private JFrame contentPanel;
	private JPanel treePanel;
	private JPanel searchButtonControlPanel;
	private JPanel resultPanel;
    private JPanel plotterPanel;
//	private JScrollPane plotterPanel; 
	private JTextField m_searchText;
	private JScrollPane jtreeFullScrollPanel;
//	private JScrollPane jtreeChoppedScrollPanel;
	private JScrollPane jtreeRecursivelyScrollPanel;
	private JScrollPane searchResultScrollPane;
	private JTextArea searchResultTextArea;
	private final int CONTENT_PANEL_WIDTH=935, CONTENT_PANEL_HEIGHT=852;
	private final int SEARCH_PANEL_WIDTH=630, SEARCH_PANEL_HEIGHT=100;
	private final int RESULT_TEXTAREA_ROWS=20, RESULT_TEXTAREA_COLUMNS=100;
	private final int RESULT_SCROLLPANEL_WIDTH=300,RESULT_SCROLLPANEL_HEIGHT=300;
	/** buttons utilities */
	private JButton mButtonFirstSearch = null;
	private JButton mButtonPrintTree = null;
	private JButton mButtonFullSearch = null;
	private JButton mButtonContinue = null;
	private JButton mButtonRefresh = null;
	private JButton	mClicked = null;
	private ActionListener mButtonListener = new ButtonListener();
	/** variables for the Search tree */
	private Enumeration<?> enume;
	/**
	 * simple constructor
	 */
	public YamlWindowPoper() {
		
	}
	/**
	 * The constructor of the class
	 * @param timeYamlFullTree
	 * @param treeNodeTimeYamlParse_RecursivelyTree_root
	 * @param timeYamlRecursivelyTree
	 * @param m_modelFullTree
	 */
	public YamlWindowPoper(JTree timeYamlFullTree, DefaultMutableTreeNode treeNodeTimeYamlParse_RecursivelyTree_root, JTree timeYamlRecursivelyTree, DefaultTreeModel m_modelRecursivelyTree) {
		this.timeYamlFullTree = timeYamlFullTree;
		this.timeYamlRecursivelyTree = timeYamlRecursivelyTree;
		this.treeNodeTimeYamlParse_RecursivelyTree_root = treeNodeTimeYamlParse_RecursivelyTree_root;
		this.m_modelRecursivelyTree = m_modelRecursivelyTree;
	}
	/**
	 * The method to pop the window panel needing the JTrees
	 */
	public void WindowPoper() {
        /**The GUI interface of the tree*/
		
		//creating the main panel
        contentPanel = new JFrame();

        //creating the tree panel
        treePanel = new JPanel();

        //creating the scroll panel to display the trees
        timeYamlFullTree.setBorder(BorderFactory.createTitledBorder("Manually generated tree"));
        jtreeFullScrollPanel = new JScrollPane(timeYamlFullTree);
        jtreeFullScrollPanel.setName("jtreeFullScrollPanel");
        jtreeFullScrollPanel.setViewportView(timeYamlFullTree);
        
/**
        timeYamlIterativelyTree.setBorder(BorderFactory.createTitledBorder("Iteratively generated tree"));
        jtreeChoppedScrollPanel = new JScrollPane(timeYamlIterativelyTree);
        jtreeChoppedScrollPanel.setName("jtreeChoppedScrollPanel");
        jtreeChoppedScrollPanel.setViewportView(timeYamlIterativelyTree);
*/
        
        timeYamlRecursivelyTree.setBorder(BorderFactory.createTitledBorder("Recursively generated tree"));
        jtreeRecursivelyScrollPanel = new JScrollPane(timeYamlRecursivelyTree);
        jtreeRecursivelyScrollPanel.setName("jtreeRecursivelyScrollPanel");
        jtreeRecursivelyScrollPanel.setViewportView(timeYamlRecursivelyTree);

        
        javax.swing.GroupLayout treePanelLayout = new javax.swing.GroupLayout(treePanel);
        treePanel.setLayout(treePanelLayout);
        treePanelLayout.setHorizontalGroup(
            treePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(treePanelLayout.createSequentialGroup()
                .addComponent(jtreeFullScrollPanel, javax.swing.GroupLayout.PREFERRED_SIZE, 393, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jtreeRecursivelyScrollPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 269, Short.MAX_VALUE)
                .addContainerGap())
        );
        treePanelLayout.setVerticalGroup(
            treePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(treePanelLayout.createSequentialGroup()
                .addGroup(treePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jtreeRecursivelyScrollPanel, 0, 0, Short.MAX_VALUE)
                    .addComponent(jtreeFullScrollPanel, javax.swing.GroupLayout.PREFERRED_SIZE, 410, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap())
        );
        
        
        treePanel.setBorder(BorderFactory.createTitledBorder(null, "Tree panel"));
        treePanel.setName("treePanel");

        //contentPanel.getContentPane().add(jtreeFullScrollPanel, BorderLayout.NORTH);
        contentPanel.getContentPane().add(treePanel, BorderLayout.NORTH);

        //creating the search panel with it buttons
        searchButtonControlPanel = new JPanel(); 
        searchButtonControlPanel.setBorder(BorderFactory.createTitledBorder(null, "Search panel")); 

        m_searchText = new JTextField(20);
        m_searchText.setToolTipText("enter the text to be searched for in the tree data structure");
        searchButtonControlPanel.add(m_searchText);

        mButtonFirstSearch = GuiUtils.makeJButton_3("Search for first occurrence", mButtonListener);
        searchButtonControlPanel.add(mButtonFirstSearch);
        mButtonPrintTree = GuiUtils.makeJButton_2("Print the tree", mButtonListener);
        searchButtonControlPanel.add(mButtonPrintTree);
        mButtonFullSearch = GuiUtils.makeJButton_2("Full search", mButtonListener);
        searchButtonControlPanel.add(mButtonFullSearch);
        mButtonRefresh = GuiUtils.makeJButton("Refresh", mButtonListener);
    	searchButtonControlPanel.add(mButtonRefresh);
        mButtonContinue = GuiUtils.makeJButton("Continue", mButtonListener);
    	searchButtonControlPanel.add(mButtonContinue);

    	searchButtonControlPanel.setSize(SEARCH_PANEL_WIDTH, SEARCH_PANEL_HEIGHT);

        //adding the search panel to the content Panel
        contentPanel.getContentPane().add(searchButtonControlPanel, BorderLayout.CENTER);

        //creating the result panel
        resultPanel = new JPanel();
        resultPanel.setBorder(BorderFactory.createTitledBorder("result Panel"));
        resultPanel.setName("resultPanel");

        
/**
        IAnalysisFactory af = IAnalysisFactory.create();
        ITree tree = af.createTreeFactory().create();
        IHistogramFactory hf = af.createHistogramFactory(tree);
        IHistogram1D h1d = hf.createHistogram1D("Test", 50, -3, 3);

        //Fill with random numbers
        Random rand = new Random();
        for (int i = 0; i < 10000; i++) h1d.fill(rand.nextGaussian());

        // Create an IPlotter
        IPlotter plotter = af.createPlotterFactory().create("Plot");
        plotter.createRegion().plot(h1d);//2.0, 2.0);
        plotter.show();
*/

        searchResultTextArea = new JTextArea();
        searchResultTextArea.setName("searchResultTextArea");

        searchResultTextArea.setRows(RESULT_TEXTAREA_ROWS);
        searchResultTextArea.setColumns(RESULT_TEXTAREA_COLUMNS);
		searchResultTextArea.setEditable(false);
		searchResultTextArea.setFont(new Font("times new roman", Font.BOLD, 11));
		
        searchResultScrollPane = new JScrollPane(searchResultTextArea);
        
		searchResultScrollPane.setName("searchResultScrollPane");
        searchResultScrollPane.setBorder(BorderFactory.createTitledBorder(null, "Search result ScrollPane"));
        searchResultScrollPane.setSize(RESULT_SCROLLPANEL_WIDTH, RESULT_SCROLLPANEL_HEIGHT);

/**        
        PieDataset dataset = createDataset();
        JFreeChart chart = createChart(dataset, "INIT classes");
        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new Dimension(339, 263));
        
        JFrame frame = new JFrame("Embedded AIDA");
        frame.getContentPane().add(chartPanel,BorderLayout.CENTER);
        frame.pack();
*/        


        plotterPanel = new JPanel();
        plotterPanel.setBorder(BorderFactory.createTitledBorder("Plotter panel"));
        plotterPanel.setName("plotterPanel");
//        plotterPanel.add(chartPanel);//,BorderLayout.CENTER);
        
        javax.swing.GroupLayout plotterPanelLayout = new javax.swing.GroupLayout(plotterPanel);
        plotterPanel.setLayout(plotterPanelLayout);
        plotterPanelLayout.setHorizontalGroup(
            plotterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 339, Short.MAX_VALUE)
        );
        plotterPanelLayout.setVerticalGroup(
            plotterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 263, Short.MAX_VALUE)
        );

        javax.swing.GroupLayout resultPanelLayout = new javax.swing.GroupLayout(resultPanel);
        resultPanel.setLayout(resultPanelLayout);
        resultPanelLayout.setHorizontalGroup(
            resultPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, resultPanelLayout.createSequentialGroup()
                .addComponent(searchResultScrollPane, javax.swing.GroupLayout.DEFAULT_SIZE, 319, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(plotterPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
        );
        resultPanelLayout.setVerticalGroup(
            resultPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(plotterPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .addComponent(searchResultScrollPane, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, 290, Short.MAX_VALUE)
        );
        
        contentPanel.getContentPane().add(resultPanel, BorderLayout.SOUTH);

        //setting the dimensions of the main panel
        contentPanel.setSize(CONTENT_PANEL_WIDTH,CONTENT_PANEL_HEIGHT);
        contentPanel.setVisible(true);  
		
	}
	/**
	 * class to attache action to actionListeners to buttons in JFrame. 
	 * @author Frederic Bonnet 
	 *
	 */
	public class ButtonListener implements ActionListener {
		public void actionPerformed(ActionEvent e) {
			if(e.getSource().equals(mButtonFirstSearch)) doFirstEnumSearchTree();
			if(e.getSource().equals(mButtonPrintTree)  ) doFullEnumPrintTree();
			if(e.getSource().equals(mButtonFullSearch) ) doFullEnumSearchTree();
			if(e.getSource().equals(mButtonRefresh)    ) doRefreshTextArea();
			if(e.getSource().equals(mButtonContinue)   ) doContinue();
		}
	}
	/**
	 * method to perform the continue action in the buttons
	 */
	protected void doContinue() {
    	mClicked = mButtonContinue;
    	contentPanel.setVisible(false);
    	contentPanel.dispose();
	}
	/**
	 * method to refresh the textArea for the mButtonRefresh
	 */
	public void doRefreshTextArea() {
		searchResultTextArea.setText(null);	
	}
	/**
	 * method to print the full tree to the textArea
	 */
	public void doFullEnumPrintTree() {
		mClicked = mButtonFullSearch;
		searchResultTextArea.append("doing full search tree for: "+m_searchText.getText()+"\n");
        System.out.println("\n");
//		enume = treeNodeTimeYamlParse_FullTree_root.preorderEnumeration();//breadthFirstEnumeration();
        enume = treeNodeTimeYamlParse_RecursivelyTree_root.preorderEnumeration();
		if ( enume != null ) {
			while ( enume.hasMoreElements() ) {
				searchResultTextArea.append(enume.nextElement().toString()+"\n");
			}
		}
	}
	/**
	 * method to perform a full search through the tree and report
	 * all of the occurrences in the tree and reports results in the 
	 * textArea panel 
	 */
	public void doFullEnumSearchTree() {
		mClicked = mButtonFullSearch;
        DefaultMutableTreeNode nodeSearchText = searchAllNode(m_searchText.getText()); 
		searchResultTextArea.append("doing full search tree for: "+m_searchText.getText()+"\n");
        DefaultMutableTreeNode node = null; 
//		enume = treeNodeTimeYamlParse_FullTree_root.preorderEnumeration();//breadthFirstEnumeration();
        enume = treeNodeTimeYamlParse_RecursivelyTree_root.preorderEnumeration();
        if (nodeSearchText != null ) {
			if ( enume != null ) {
				while ( enume.hasMoreElements() ) {
					node = (DefaultMutableTreeNode)enume.nextElement();
					if (m_searchText.getText().equals(node.getUserObject().toString())) {
						//ArrayList<Double> numbers = new ArrayList();
						//numbers.add(Double.parseDouble(enume.nextElement().toString()));
						searchResultTextArea.append(enume.nextElement().toString());//+" elements of list"+enume.nextElement()+"\n");
			            //make the node visible by scroll to it 
						TreeNode[] nodes = m_modelRecursivelyTree.getPathToRoot(node);
				        TreePath path = new TreePath(nodes); 		
						timeYamlFullTree.scrollPathToVisible(path);  
			            timeYamlFullTree.setSelectionPath(path);
					}
				}
			}
		}else {
			//node with string not found show message 
			JOptionPane.showMessageDialog(this,  
	        "Node with string " + m_searchText.getText() + " not found",  
	        "Node not found", JOptionPane.INFORMATION_MESSAGE);                          
		}
	}
	/**
	 * method to search for the existence of the string in the tree 
	 * @param nodeStr
	 * @return
	 */
	public DefaultMutableTreeNode searchAllNode(String nodeStr) {
		DefaultMutableTreeNode node = null;
        //Get the enumeration
        //enume = treeNodeTimeYamlParse_FullTree_root.preorderEnumeration();
		enume = treeNodeTimeYamlParse_RecursivelyTree_root.preorderEnumeration();
        //iterate through the enumeration 
        while(enume.hasMoreElements()) 
        { 
            //get the node
            node = (DefaultMutableTreeNode)enume.nextElement();             
            //match the string with the user-object of the node 
            if(nodeStr.equals(node.getUserObject().toString())) 
            {
                return node;
            }
        }         
        //tree node with string node found return null 
        return null; 
	}
	/**
	 * method to perform a search for the first occurrence within the 
	 * tree and report to the textArea the result. 
	 */
	public void doFirstEnumSearchTree() {
		mClicked = mButtonFirstSearch;
        DefaultMutableTreeNode node = searchNode(m_searchText.getText()); 
        
        if(node != null) 
        { 
            //make the node visible by scroll to it 
            TreeNode[] nodes = m_modelRecursivelyTree.getPathToRoot(node);
            TreePath path = new TreePath(nodes); 
            timeYamlFullTree.scrollPathToVisible(path);  
            timeYamlFullTree.setSelectionPath(path);
            timeYamlFullTree.makeVisible(path);
        } 
        else 
        { 
            //node with string not found show message 
            JOptionPane.showMessageDialog(this,  
            "Node with string " + m_searchText.getText() + " not found",  
            "Node not found", JOptionPane.INFORMATION_MESSAGE);                          
        } 
	}
    /** 
     * This method takes the node string and 
     * traverses the tree till it finds the node 
     * matching the string. If the match is found  
     * the node is returned else null is returned 
     *  
     * @param nodeStr node string to search for 
     * @return tree node  
     */ 
    public DefaultMutableTreeNode searchNode(String nodeStr) 
    {
        DefaultMutableTreeNode node = null;
        //Get the enumeration
        //enume = treeNodeTimeYamlParse_FullTree_root.breadthFirstEnumeration();
        enume = treeNodeTimeYamlParse_RecursivelyTree_root.breadthFirstEnumeration();
        //iterate through the enumeration 
        while(enume.hasMoreElements()) 
        {
            //get the node 
            node = (DefaultMutableTreeNode)enume.nextElement(); 

            //match the string with the user-object of the node 
            if(nodeStr.equals(node.getUserObject().toString())) 
            { 
                //tree node with string found
            	searchResultTextArea.append("found a result\n");
            	searchResultTextArea.append(node.getUserObject().toString()+"\n");
            	searchResultTextArea.append(node.getLastLeaf().toString()+"\n");
                return node;                          
            }
        }    
        //tree node with string node found return null 
        return null; 
    }

}
