import React, { Component } from 'react';
import ReactTable from "react-table";
const d3 = require("d3");
import 'react-table/react-table.css'

class Overview extends Component {
  constructor(props) {
    super(props);
    this.state = { data: null };
  }
  componentDidMount() {
    d3.csv('data/cluster_information.csv', (error, data) => {
      data.forEach(row=>{
        row.visit1 = +row.visit1;
        row.visit2 = +row.visit2;
        row.visit3 = +row.visit3;
      });
      this.setState({data: data});
    });
  }
  render() {
    return (<div>
      {this.state.data ? <ReactTable
        data={this.state.data}
        columns={[
          { Header: "Patient ID", accessor: 'patient_id'},
          { Header: "V-gene allele", accessor: 'vgene_allele'},
          { Header: "Visit 1 size", accessor: 'visit1'},
          { Header: "Visit 2 size", accessor: 'visit2'},
          { Header: "Visit 3 size", accessor: 'visit3'}
        ]}
      /> : null }
    </div>);
  }
}

module.exports = Overview;
