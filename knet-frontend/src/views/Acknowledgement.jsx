import React from 'react';

const Acknowledgement = () => (
    <div className="container">
    <div className="text-section">
      <h3>We would like to thank the following individuals and organizations for their invaluable contributions to this project:</h3>
      <ul>
          <li>Project Contributors: Jiaming Hu, Lukas Rupprecht, and Stephan Schürer</li>
          <br></br>
          <li>Funding and Support:</li>
          <ul>
            <li> This work was supported by National Institutes of Health (NIH) grants U24TR002278 (Illuminating the Druggable Genome) administered by NCATS, the data curation projects U01LM012630 and R01LM01339 from the National Library of Medicine (NLM), and by the State of Florida Bankhead-Coley Cancer Research Program, grant 23B16.</li>
          </ul>
          <br></br>
          <li>Tools and Resource:
            <ul>
                <li>ChEMBL, a comprehensive chemical database of drug-like molecules, https://www.ebi.ac.uk/chembldb/</li>
                <li>KKB, a large-scale database of kinase structure-activity data, http://www.eidogen.com/kinasekb.php</li>
                <li>Rdkit, an open-source toolkit for cheminformatics</li>
            </ul>
          </li>
      </ul>
      <br></br>
      <br></br>
      <h3>Please acknowledge KNET in your publications by citing the following reference:</h3>
      <h3>Jiaming Hu, Lukas Rupprecht, and Stephan C. Schürer:  Multitask Deep Neural Network for Kinase Classification (KNET) 2024 https://knet.idsc.miami.edu/</h3>
      <br></br>
      <h3>If you have any questions, feel free to reach out to us at sschurer@miami.edu</h3>
    </div>
</div>
);

export default Acknowledgement;