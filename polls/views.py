from django.shortcuts import get_object_or_404, render
from django.http import HttpResponse
from .models import Question
from django.http import Http404
from .forms import NameForm
from django.http import HttpResponseRedirect


def index(request):
    latest_question_list = Question.objects.order_by('-pub_date')[:5]
    context = {'latest_question_list': latest_question_list}
    return render(request, 'polls/index.html', context)


def detail(request, question_id):
    question = get_object_or_404(Question, pk=question_id)
    return render(request, 'polls/detail.html', {'question': question})


def results(request, question_id):
    response = "You're looking at the results of question %s."
    return HttpResponse(response % question_id)


def vote(request, question_id):
    return HttpResponse("You're voting on question %s." % question_id)


def gcam(request):
    # if this is a POST request we need to process the form data
	if request.method == 'POST':
		# create a form instance and populate it with data from the request:
		form = NameForm(request.POST)
		# check whether it's valid:
		if form.is_valid():
			name = request.POST.get('your_name','')
			print('Entered name:', name)
			return HttpResponseRedirect('redirect/')

    # if a GET (or any other method) we'll create a blank form
	else:
		form = NameForm()

	return render(request, 'polls/gcam.html', {'form': form})

def redirect(request):
    return render(request, 'polls/redirect.html')
