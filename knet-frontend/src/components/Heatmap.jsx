import React from 'react';
import Plot from 'react-plotly.js';

const Heatmap = ({ data }) => {
  if (data.length === 0) {
    return <p>No data provided</p>;
  }

  // Extract `standardizedSmiles` for y-axis labels
  const yAxisLabels = data.map((item) => item.standardizedSmiles);

  // Extract columns, excluding non-relevant keys
  const columns = Object.keys(data[0]).filter(
    (key) => key !== 'standardizedSmiles' && key !== 'originalSmiles' && key !== 'personalId'
  );

  // Construct the data matrix for the heatmap
  const matrix = data.map((item) => columns.map((col) => item[col]));

  return (
    <div style={{ width: '100%' }}>
      <Plot
        data={[
          {
            z: matrix,
            x: columns,
            y: yAxisLabels,
            type: 'heatmap',
            colorscale: [
              [0, 'blue'], // 0 -> blue
              [0.5, 'white'], // 0.5 -> white
              [1, 'red'], // 1 -> red
            ],
            colorbar: {
              tickvals: [0, 0.5, 1],
              ticktext: ['0', '0.5', '1'],
              x: 0.95, // Position the colorbar within the visible area
              xanchor: 'left',
              lenmode: 'fraction',
              len: 1.0,
            },
          },
        ]}
        layout={{
          title: 'Heatmap',
          height: 600, // Adjust the height to your preference
          autosize: true, // Let Plotly handle the sizing
          xaxis: {
            domain: [0, 0.9], // Adjusted domain to leave space for the colorbar
            tickmode: 'array',
            tickvals: columns,
            tickangle: 45, // Rotate x-axis labels for better readability
            automargin: true,
          },
          yaxis: {
            automargin: true, // Automatically adjust margins to prevent cutoff
          },
          margin: {
            l: 150, // Increase left margin to fit longer y-axis labels
            r: 50,
            b: 150, // Bottom margin for x-axis labels
            t: 50,
          },
        }}
        useResizeHandler={true}
        style={{ width: '100%', height: '600px' }}
      />
    </div>
  );
};

export default Heatmap;
