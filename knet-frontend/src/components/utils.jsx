import * as XLSX from 'xlsx';
import Papa from 'papaparse';

export const exportToCSV = (dataToExport, fileName) => {
  // Check if there's at least one non-empty 'personalId'
  const hasPersonalId = dataToExport.some(row => row.personalId && row.personalId.trim() !== '');

  // Remove 'personalId' from each row if it's empty and prepare for export
  const exportData = dataToExport.map(row => {
    const { personalId, ...rest } = row;
    if (hasPersonalId) {
      // Move personalId to the first position if present
      return { personalId, ...rest };
    }
    // Exclude personalId if all are empty
    return rest;
  });

  // Convert data to a worksheet
  const worksheet = XLSX.utils.json_to_sheet(exportData);

  // Convert worksheet to CSV format
  const csvContent = XLSX.utils.sheet_to_csv(worksheet);

  // Create a blob and trigger download
  const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
  const link = document.createElement('a');
  const url = URL.createObjectURL(blob);

  link.setAttribute('href', url);
  link.setAttribute('download', fileName || 'export.csv');
  link.style.visibility = 'hidden';

  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
};

export const parseCsvFile = (event) => {
  return new Promise((resolve, reject) => {
    const file = event.target.files[0];
    if (file) {
      Papa.parse(file, {
        header: false,
        skipEmptyLines: true, // Skip empty lines
        complete: (results) => {
          const data = results.data;
          const mappedData = data.map(row => ({
            personal_id: Object.values(row)[0]?.trim(), // First column
            smiles_code: Object.values(row)[1]?.trim(), // Second column
          }));
          resolve(mappedData);
        },
        error: (error) => {
          console.error('Error parsing CSV:', error);
          reject(error);
        },
      });
    } else {
      reject(new Error('No file provided'));
    }
  });
};
