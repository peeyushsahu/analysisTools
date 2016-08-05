
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