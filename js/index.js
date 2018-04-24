import React, { Component } from 'react';
import ReactDOM from 'react-dom';
import { BrowserRouter as Router, Route } from 'react-router-dom';
import { LinkContainer } from 'react-router-bootstrap';
import { Nav, Navbar, NavItem, Grid } from 'react-bootstrap';
import 'bootstrap/dist/css/bootstrap.css';

import Trees from './trees.js';
import Alignments from './alignments.js';
import JSONs from './jsons.js';


class App extends Component {
  render(){
    return(<Router>
      <div>
        <Navbar>
          <Navbar.Header>
            <Navbar.Brand>
              VEG Immunology Visualization
            </Navbar.Brand>
          </Navbar.Header>
          <Nav>
            <LinkContainer exact to="/">
              <NavItem eventKey={1} href="#">
                Trees
              </NavItem>
            </LinkContainer>
            <LinkContainer to="/alignments">
              <NavItem eventKey={2} href="#">
                Alignments
              </NavItem>
            </LinkContainer>
            <LinkContainer to="/jsons">
              <NavItem>
                JSONs
              </NavItem>
            </LinkContainer>
          </Nav>
        </Navbar>
        <Grid>
          <Route exact path="/" component={Trees} />
          <Route path="/alignments" component={Alignments} />
          <Route path="/jsons" component={JSONs} />
        </Grid>
      </div>
    </Router>);
  }
}

ReactDOM.render(
  <App />,
  document.body.appendChild(document.createElement('div'))
)
