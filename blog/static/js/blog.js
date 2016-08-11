
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

function readtable4txt(file_path){
    d3.text(file_path, function(error, data) {
    console.log('This is here');
    console.log(d3.tsv.parseRows(data));
    });
}

$(document).ready(function() {
    $('#example').DataTable( {
        columns: [
            { title: "celltype" },
            { title: "genecluster" },
            { title: "p-val" },
            { title: "FDR" },
        ]
    } );
} );