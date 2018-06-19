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
    this.state = { active: "no_replicates", size: 200 };
  }
  onSelect(key){
    if(key) {
      const active = key.split('-')[0],
        size = +key.split('-')[1];
      this.setState({
        active: active,
        size: size
      });
    }
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
              <MenuItem eventKey='replicates-30' active={this.state.active == 'replicates' && this.state.size == 30}>
                Replicates, size 30
              </MenuItem>
              <MenuItem eventKey='no_replicates-30' active={this.state.active == 'no_replicates' && this.state.size == 30}>
                No replicates, size 30
              </MenuItem>
              <MenuItem eventKey='replicates-200' active={this.state.active == 'replicates' && this.state.size == 200}>
                Replicates, size 200
              </MenuItem>
              <MenuItem eventKey='no_replicates-200' active={this.state.active == 'no_replicates' && this.state.size == 200}>
                No replicates, size 200
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
          <Route exact path="/" render={ () => <Trees active={this.state.active} size={this.state.size} /> } />
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
