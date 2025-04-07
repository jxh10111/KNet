import React from 'react';
import { BrowserRouter as Router, Route, Routes, Link } from 'react-router-dom';
import About from './views/About';
import Standardizer from './views/Standardizer';
import Submission from './views/Submission';
import Acknowledgement from './views/Acknowledgement';
import { AppBar, Toolbar, Typography, Button, Container } from '@mui/material';

import './App.css';

const App = () => {
  return (
<Router>
      <AppBar position="fixed">
        <Toolbar>
          <Typography variant="h6" sx={{ flexGrow: 1 }}>
            KNET
          </Typography>
          <Button color="inherit" component={Link} to="/">Home</Button>
          <Button color="inherit" component={Link} to="/submit">Predictor</Button>
          <Button color="inherit" component={Link} to="/standardizer">Standardizer</Button>
          <Button color="inherit" component={Link} to="/acknowledgement">Acknowledgement</Button>
        </Toolbar>
      </AppBar>

      <Container> {/* Add margin to avoid overlap with AppBar */}
        <Routes>
          <Route path="/" element={<About />} />
          <Route path="/submit" element={<Submission />} />
          <Route path="/standardizer" element={<Standardizer />} />
          <Route path="/acknowledgement" element={<Acknowledgement />} />
        </Routes>
      </Container>
    </Router>
  );
};

export default App;
