import React from 'react';
const About = () => (
      <div className="container">
      <img
        src="bitmap-new.png"
        alt="Top Image"
        className="top-image"
      />
      <div className="text-section">
        <h3>Welcome to our compound-kinase activity prediction platform! This application allows you to predict the probability of inhibition for small molecules, represented by SMILES strings, against a broad range of protein kinases.</h3>
        <h3>
        The predictions are powered by two multitask deep neural network (MTDNN) models trained using large-scale bioactivity data with over 1 million human kinase bioactivity annotations and more than 400,000 unique small molecules from databases ChEMBL and the Kinase Knowledge Base (KKB): </h3>
        <ul>
            <li>Kinase-Cutoff6: This model covers 406 kinases and classifies molecules as active or inactive based on a pActivity cutoff value of 6. It provides wide kinase coverage, offering broad insights across many targets.</li>
            <br></br>
            <li>Kinase-Cutoff7: This model covers 328 kinases and classifies molecules using a stricter pActivity cutoff value of 7. While the kinase coverage is smaller, it provides higher accuracy and more reliable predictions. It is ideal for users seeking precise inhibition probabilities for a refined set of kinases.</li>
        </ul>
      </div>
  </div>
);

export default About;