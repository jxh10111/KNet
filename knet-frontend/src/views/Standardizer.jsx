import React, { useState } from 'react';
import { TextField, Button, Box, Typography } from '@mui/material';
import axios from 'axios';
import Grid from '@mui/material/Grid';
import Input from '@mui/material/Input';
import {exportToCSV, parseCsvFile } from '../components/utils'

const Standardizer = () => {
  const [smilesInput, setSmilesInput] = useState('');
  const [validationResults, setValidationResults] = useState([]);
  const [uploadedData, setUploadedData] = useState([]);
  const [isLoading, setIsLoading] = useState(false);
  const [errorMessage, setErrorMessage] = useState('');



  const handleValidation = async () => {
    setIsLoading(true);
    const updatedMapping = uploadedData;
    console.log(updatedMapping)
    const inputSmiles = [];
    const inputLines = smilesInput.split('\n').map(s => s.trim()).filter(s => s);
    console.log(inputLines)
    inputLines.forEach(line => {
      const parts = line.split(';');
      let smilesCode = '';
      let personalId = '';
  
      if (parts.length === 2) {
        personalId = parts[0].trim();
        smilesCode = parts[1].trim();
      } else {
        console.warn(`Invalid line format: ${line}`);
        return;
      }
      if (smilesCode) {
        inputSmiles.push(smilesCode);
        updatedMapping.push({ smiles_code:smilesCode, personal_id:personalId }); // Map SMILES code to personalId
      }
    });
    
    try {
      const response = await axios.post('/standardize-smiles', {
        smiles: updatedMapping
      });
      setValidationResults(response.data);
    } catch (error) {
      console.error('Error validating SMILES:', error);
      setErrorMessage(error.message);

    } finally {
      setIsLoading(false);
    };
  };

  const handleFileUpload = async (event) => {
    const mappedData = await parseCsvFile(event);
    setUploadedData(mappedData);
    };
  


  return (
    <Box className="container" sx={{ mt: 2 }}>
        <div className="text-section">
        <h3>
        This tool allows you to standardize SMILES strings using a standardized pipeline to ensure they are correctly formatted. 
        You can submit SMILES strings in one of the following ways:
        </h3>
        <ol>
          <li>
            Copy & Paste: Copy and paste the SMILES strings along with a required corresponding identifiers into the text box.
            Each SMILES string should be on a new line, with the identifier and SMILES separated by a semicolon (;).
            <br />
            <strong>Example:</strong> <code>Molecule1;CSc1cn(C2O[C@H](C)[C@@H](O)[C@H]2O)c3ncnc(N)c13</code>
          </li>
          <li>
          Upload a .csv file:
          <br />
          <strong>Example:</strong>
          <img src="/csvContent.png" alt="Example CSV content" width="100%"></img>
          </li>
        </ol>
      </div>
      <Box sx={{ mt: 2 }}>
        <TextField
          label="Enter SMILES strings, one per line (REQUIRED with personalized ID; SMILES)"
          placeholder="Molecule1;CSc1cn(C2O[C@H](C)[C@@H](O)[C@H]2O)c3ncnc(N)c13"
          variant="outlined"
          fullWidth
          multiline
          rows={4}
          value={smilesInput}
          onChange={(e) => setSmilesInput(e.target.value)}
        />
      </Box>
      <Grid item xs={12} sm={3}>
          <Input
            type="file"
            inputProps={{ accept: '.csv' }}
            onChange={handleFileUpload}
          />
        </Grid>

      <Box sx={{ mt: 2 }}>
        <Button variant="contained" onClick={handleValidation}>
          Standardize SMILES
        </Button>
      </Box>

      {uploadedData.length > 0 && (
        <div>
          <p>{uploadedData.length} records loaded from CSV file.</p>
        </div>
      )}
            {isLoading && (
        <div className="loading-indicator">
          <div className="spinner"></div>
        </div>
      )}
            {errorMessage && <div className="error-message">An unexpected error occurred. Validate that the SMILES codes are valid and try again.</div>}

      {validationResults.length > 0 && (
      <div>
      <Grid container spacing={2} direction="column" className="table-section">
      <Grid item>
      <Typography variant="h6">Standardized Data:</Typography>
      </Grid>
      <Grid item>
      <Button variant="contained" color="primary" onClick={() => exportToCSV(validationResults, "standardized_smiles.csv")}>
        Export as CSV
      </Button>
      </Grid>
      </Grid>
      </div>

)}

      {validationResults.length > 0 && (
        <Box sx={{ mt: 2 }}>
        <Typography variant="h6">Validation Results:</Typography>
        {validationResults.map(({ originalSmiles, isValid, ...attributes }, index) => (
          <Box key={index} sx={{ mb: 1 }}>
            <Typography sx={{ color: isValid ? 'green' : 'red' }}>
              {`SMILES: ${originalSmiles} - ${isValid ? 'Valid' : `Invalid`}`}
            </Typography>
            {isValid && (
              <Box sx={{ ml: 2 }}>
                {Object.entries(attributes).map(([key, value]) => (
                  <Typography key={key} sx={{ color: 'green' }}>
                    {`${key}: ${value !== null ? value : 'N/A'}`}
                  </Typography>
                ))}
              </Box>
            )}
          </Box>
        ))}
      </Box>
      )}
    </Box>
  );
};

export default Standardizer;
