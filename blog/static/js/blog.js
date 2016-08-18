
function sayBye(){
    $('#title').html('Bye');
    $('#title').click(function(){
        $('#title').html('Hello');
        $('#title').off('click');
    });
};


function dim(bool){
    if (typeof bool=='undefined') bool=true; // so you can shorten dim(true) to dim()
    document.getElementById('dimmer').style.display=(bool?'block':'none');
}

function readtable(file_path, divid){
    d3.text(file_path, function(error, data){
        // Creates the DOM elements table, tbody and thead
        var TABLE = d3.select(divid).append("table")
                .attr("class","display")
                .attr("width","100%")
                .attr("id","example")
                .style("visibility", "hidden"),

            thead = TABLE.append("thead")
            tbody = TABLE.append("tbody");

        var parsedData = d3.tsv.parseRows(data);
        // .shift() removes the first element of the array parsedData and modifies the array.
        // It returns this element and I store it into the variable called header
        var header = parsedData.shift();

        // Create the rows to the table based on the number of entries in the parsedData array
        var rows = tbody.selectAll("tr")
            .data(parsedData, String)
            .enter()
            .append("tr")

        // Append a child row element "tr" to <thead> where to put the header
        var rowsH = thead.append("tr");

        // Create the cells with <td> for each element of the array of array parsedData
        var cells = rows.selectAll("td")
            .data(function(c, i) {
                return parsedData[i];
            })
            .enter()
            .append("td")
            .text(function(d) { return d; })
            // Adding a class element corresponding to the array position for column highlight on mouseover
            .attr('class', function(d, i){ return "col_" + i; });

        // Creates the header cell
        var cellsH = rowsH.selectAll("th")
            .data(header)
            .enter()
            .append("th")
            .text(function(d) { return d; })
            .attr('class', function(d, i){ return "colH_" + i; });

        // Make the table zebra coloured
        tbody.selectAll("tr").style("background-color", function(d, i){
            if(i%2 == 0){
                return "#D8D8D8";
            }
        })

    $(document).ready(function() {
        $('#example').DataTable( {
            columns: [
                { title: "celltype" },
                { title: "genecluster" },
                { title: "p-val" },
                { title: "FDR" }
            ]
        } );
    } );
    TABLE.style("visibility", "visible");
    });
}


function ifsynonym(){
    //console.log(document.getElementById("organism").style.visibility);
    var state = document.getElementById("organism").style.visibility
    if (state == 'visible')
        document.getElementById("organism").style.visibility="hidden";
    if (state == 'hidden')
        document.getElementById("organism").style.visibility="visible";
}


function show_control_sample_name(){
    var theSelect = ExprbasedForm.select_control;
    var value = theSelect[theSelect.selectedIndex].value;
    //console.log(theSelect[theSelect.selectedIndex].value);
    if (value == 'sample')
        $('#id_control_sample_name').attr('disabled', false);
    if (value == 'mean')
        $('#id_control_sample_name').attr('disabled', true);
}


$(window).scroll(function() {
    if( $(this).scrollTop() > 150 ) {
        $(".navbar").addClass("navbar-scrolled");
    } else {
        $(".navbar").removeClass("navbar-scrolled");
    }
});


function upload_example_genebased(){
    $('#id_gene_text_field').val("oct4\nsox2\ntrim71\nklf4\nc-myc\nn-myc\ncdkn1a\nAATF\nski\nRUNX1\np300");
}








