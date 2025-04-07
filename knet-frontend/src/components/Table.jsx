import React from 'react';
import {
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow,
  Paper,
} from '@mui/material';

export const DataTable = ({ data }) => {
  if (!data || data.length === 0) {
    return <div>No data available</div>; // Handle empty data case
  }

  const columns = Object.keys(data[0]); // Extract keys from the first object as column headers

  return (
    <TableContainer component={Paper}>
      <Table>
        <TableHead>
          <TableRow>
            {columns.map((key, index) => (
              <TableCell key={index}>{key}</TableCell>
            ))}
          </TableRow>
        </TableHead>
        <TableBody>
          {data.map((row, rowIndex) => (
            <TableRow key={rowIndex}>
              {columns.map((column, colIndex) => (
                <TableCell key={colIndex}>
                  {typeof row[column] === 'object'
                    ? JSON.stringify(row[column])
                    : row[column]}
                </TableCell>
              ))}
            </TableRow>
          ))}
        </TableBody>
      </Table>
    </TableContainer>
  );
};
