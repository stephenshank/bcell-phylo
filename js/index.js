import React, { Component } from 'react';
import ReactDOM from 'react-dom';
import { BrowserRouter as Router, Route } from 'react-router-dom';
import { LinkContainer } from 'react-router-bootstrap';
import { Nav, Navbar, NavItem, NavDropdown, MenuItem, Grid } from 'react-bootstrap';
import 'bootstrap/dist/css/bootstrap.css';

import Trees from './trees.js';
import JSONs from './jsons.js';

require('./main.css');

class App extends Component {
  constructor(props){
    super(props);
    this.state = { active: "no_replicates" };
  }
  onSelect(key){
    if(key) this.setState({ active: key });
  }
  render(){
    return(<Router>
      <div>
        <Navbar onSelect={key=>this.onSelect(key)} fixedTop>
          <Navbar.Header>
            <Navbar.Brand>
              VEG Immunology Visualization
            </Navbar.Brand>
          </Navbar.Header>
          <Nav>
            <NavDropdown title='Dataset'>
              <MenuItem eventKey='replicates' active={this.state.active == 'replicates'}>
                Replicates
              </MenuItem>
              <MenuItem eventKey='no_replicates' active={this.state.active == 'no_replicates'}>
                No replicates
              </MenuItem>
            </NavDropdown>
            <LinkContainer exact to="/">
              <NavItem href="#">
                Trees
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
          <Route exact path="/" render={()=><Trees active={this.state.active} />} />
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
