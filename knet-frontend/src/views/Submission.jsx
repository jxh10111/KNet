import React, { useState } from 'react';
import Button from '@mui/material/Button';
import Input from '@mui/material/Input';

import InputLabel from '@mui/material/InputLabel';
import MenuItem from '@mui/material/MenuItem';
import FormControl from '@mui/material/FormControl';
import Select from '@mui/material/Select';
import { Typography, Box } from '@mui/material';
import TextField from '@mui/material/TextField';
import Grid from '@mui/material/Grid';
import Heatmap from '../components/Heatmap';

import { exportToCSV, parseCsvFile } from '../components/utils';

import { DataTable } from '../components/Table';

const Submission = () => {
  const [inputText, setInputText] = useState('');
  const [selectedModel, setSelectedModel] = useState('kinase-cutoff7');
  const [isLoading, setIsLoading] = useState(false);
  const [errorMessage, setErrorMessage] = useState('');
  const [dataMatrix, setDataMatrix] = useState([]);
  const [uploadedData, setUploadedData] = useState([]);

  const handleFileUpload = async (event) => {
    const mappedData = await parseCsvFile(event);
    setUploadedData(mappedData);
    };

  const handleSubmit = () => {
    setIsLoading(true);
    setErrorMessage('');
    setDataMatrix([]);
  
    const updatedMapping = uploadedData;
    const inputSmiles = [];
    const inputLines = inputText.split('\n').filter(line => line.trim() !== '');
  
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
  

    // **Ensure that only SMILES codes are sent to the backend**
    const requestBody = { "smiles": updatedMapping };
    const requestHeaders = new Headers({ "Content-Type": "application/json" });
  
    // Update the URL based on the selected model
    const cutoff = selectedModel.replace('kinase-', '');
    const url = `/predict/${cutoff}`;
  
    fetch(url, {
      method: "POST",
      body: JSON.stringify(requestBody),
      headers: requestHeaders
    })
    .then((response) => {
      if (response.ok) {
        return response.json();
      } else {
        console.log([...response.headers.entries()]);
        console.log(response.body);
        throw new Error("Request failed with status code " + response.status);
      }
    })
    .then((json) => {
      const data = json.prediction;
  
      // Process unique SMILES and Proteins
  

  
      setDataMatrix(data);
      setErrorMessage('');
    })
    .catch((error) => {
      console.error(error);
      setErrorMessage(error.message);
    })
    .finally(() => {
      setIsLoading(false);
    });
  };  

  return (
    <Box className="container" sx={{ mt: 2 }}>
      <h3>To run the prediction:</h3>
      <ul>
        <li>Ensure your SMILES strings are valid using a SMILES Validator</li>
        <br />
        <li>
          Upload a <code>.csv</code> file, or paste each SMILES string on a new line.
          When submitting SMILES strings, a personal ID is required. Please use the following format:
          <br />
          <strong>Example:</strong> <code>Molecule1;CSc1cn(C2O[C@H](C)[C@@H](O)[C@H]2O)c3ncnc(N)c13</code>
        </li>
        <br />
        <li>Results will be generated and available for download once finished.</li>
      </ul>
      <Grid container spacing={2} alignItems="center">
        {/* Select Component */}
        <Grid item xs={12} sm={3}>
          <FormControl fullWidth>
            <InputLabel id="model-select-label">Model</InputLabel>
            <Select
              labelId="model-select-label"
              id="model-select"
              value={selectedModel}
              label="Selected Model"
              onChange={(e) => setSelectedModel(e.target.value)}
              sx={{ minWidth: 200 }}
            >
              <MenuItem value={'kinase-cutoff6'}>Kinase-Cutoff6</MenuItem>
              <MenuItem value={'kinase-cutoff7'}>Kinase-Cutoff7</MenuItem>
            </Select>
          </FormControl>
        </Grid>

        {/* TextField Component */}
        <Grid item xs={12} sm={6}>
          <TextField
            fullWidth
            multiline
            rows={4}
            placeholder="Molecule1;CSc1cn(C2O[C@H](C)[C@@H](O)[C@H]2O)c3ncnc(N)c13"
            id="outlined-basic"
            label="Enter SMILES strings, one per line (REQUIRED with personalized ID;SMILES Code)"
            variant="outlined"
            value={inputText}
            onChange={(e) => setInputText(e.target.value)}
            sx={{ width: '100%' }}
          />
        </Grid>

        {/* File Upload Component */}
        <Grid item xs={12} sm={3}>
          <Input
            type="file"
            inputProps={{ accept: '.csv' }}
            onChange={handleFileUpload}
          />
        </Grid>

        {/* Submit Button */}
        <Grid item xs={12}>
          <Button
            onClick={handleSubmit}
            variant="contained"
            fullWidth
            sx={{ height: '100%', mt: 2 }}
          >
            Submit
          </Button>
        </Grid>
      </Grid>

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

      {dataMatrix.length > 0 && !isLoading && (
      <div>
      <Grid container spacing={2} direction="column" className="table-section">
      <Grid item>
      <Typography variant="h6">Fetched Data:</Typography>
      </Grid>
      <Button variant="contained" color="primary" onClick={() => exportToCSV(dataMatrix, "prediction_smiles.csv")}>
        Export as CSV
      </Button>
      <DataTable data={dataMatrix} />
      <Grid item>
      </Grid>
      </Grid>
      </div>

)}
       <div>
      {dataMatrix.length > 0 && (
        <Heatmap data={dataMatrix} />
      )}
    </div>
    </Box>
  );
};

export default Submission;
