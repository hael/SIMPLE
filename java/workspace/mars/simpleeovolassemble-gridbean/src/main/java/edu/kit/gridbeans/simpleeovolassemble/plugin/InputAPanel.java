package edu.kit.gridbeans.simpleeovolassemble.plugin;

import java.util.logging.Level;
import java.util.logging.Logger;
import com.intel.gpe.clients.api.Client;
import com.intel.gpe.gridbeans.plugins.DataSetException;
import com.intel.gpe.gridbeans.plugins.swing.panels.GridBeanPanel;
import com.intel.gpe.gridbeans.plugins.translators.GPELocalFileValueTranslator;
import com.intel.gpe.gridbeans.plugins.validators.FileIsReadableAndNotDirectoryValidator;
import com.intel.gpe.util.defaults.preferences.INode;
import com.intel.gpe.util.swing.controls.configurable.TextEditor;
import edu.kit.gridbeans.simpleeovolassemble.SimpleEoVolAssembleGridBean;


public class InputAPanel extends GridBeanPanel 
{

  public InputAPanel(Client client, INode node) 
  {
    super(client, "InputA editor", node);
    initComponents();

    try 
    {
      linkComponents();
    } 

    catch (DataSetException ex) 
    {
      Logger.getLogger(InputAPanel.class.getName()).log(Level.SEVERE, null, ex);
    }
  }


  /** Links this panel to the Model */
  private void linkComponents() throws DataSetException 
  {
    TextEditor editor1 = new TextEditor(getNode(), getClient());
    add(editor1);
    linkTextEditor(SimpleEoVolAssembleGridBean.INPUT_FILE, editor1);
    setValueTranslator(SimpleEoVolAssembleGridBean.INPUT_FILE, new GPELocalFileValueTranslator());
    setValueValidator(SimpleEoVolAssembleGridBean.INPUT_FILE, FileIsReadableAndNotDirectoryValidator.getInstance());
  }

  @SuppressWarnings("unchecked")
  private void initComponents() 
  {
    setLayout(new java.awt.BorderLayout());
  }

}
